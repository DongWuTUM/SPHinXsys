/**
* @file 	dw_3d_pinched_hemisphere.cpp
* @brief 	This is the benchmark test of the shell.
* @details  We consider large deformation of a hemisphere surface.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
* @version  0.1
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real radius = 9.98;                                     /** Radius of the inner wall of the hemisphere. */
Real thickness = 0.04;                                  /** Thickness of the hemisphere. */
Real hole_angle = 18.0 / 180.0 * Pi;                    /** Angle of hole located at the hemisphere top. */
Real radius_mid_surface = radius + thickness / 2.0;      /** Radius of the mid surface. */
int particle_number = 160;								 /** Particle number in the bottom edge. */
Real particle_spacing_ref = 2.0 * Pi * radius_mid_surface / (Real)particle_number;
int particle_number_height = radius_mid_surface * (0.5 * Pi - hole_angle) / particle_spacing_ref;
int BWD = 0;
Real BW = particle_spacing_ref * (Real)BWD;              /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-radius - thickness, -radius - thickness, 0.0),
	Vec3d(radius + thickness, radius + thickness, radius + thickness));

/** For material properties of the solid. */
Real rho0_s = 1.0; 			                             /** Normalized density. */
Real Youngs_modulus = 6.825e7;	                         /** Normalized Youngs Modulus. */
Real poisson = 0.3; 			                         /** Poisson ratio. */
Real physical_viscosity = 60.0;                          /** physical damping, here we choose the same value as numerical viscosity. */

Real F_full = 400;
std::vector<Vecd> point_forces{ Vec3d(F_full, 0.0, 0.0), Vec3d(-F_full, 0.0, 0.0), 
							    Vec3d(0.0, F_full, 0.0), Vec3d(0.0, -F_full, 0.0) };
std::vector<Vecd> reference_positions{ Vec3d(radius_mid_surface, 0.0, 0.0), Vec3d(-radius_mid_surface, 0.0, 0.0), 
									   Vec3d(0.0, -radius_mid_surface, 0.0), Vec3d(0.0, radius_mid_surface, 0.0) };
Real time_to_full_external_force = 0.1;

// Observer location
StdVec<Vecd> observation_location = { Vecd(radius_mid_surface, 0.0, 0.0), 
									  Vecd(0.0, radius_mid_surface, 0.0) };

/** Define application dependent particle generator for thin structure. */
class HemisphereParticleGenerator : public SurfaceParticleGenerator
{
public:
	HemisphereParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body) {};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number_height + 1; i++)
		{
			Real circle_radius = radius_mid_surface * cos(i * (0.5 * Pi - hole_angle) / particle_number_height);
			int circle_particle_number = circle_radius * 2.0 * Pi / particle_spacing_ref;
			for (int j = 0; j < circle_particle_number; j++)
			{
				Real x = circle_radius * sin(j * 2 * Pi / (Real)circle_particle_number);
				Real y = circle_radius * cos(j * 2 * Pi / (Real)circle_particle_number);
				Real z = radius_mid_surface * sin(i * (0.5 * Pi - hole_angle) / particle_number_height);
				initializePositionAndVolumetricMeasure(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref);
				Vecd n_0 = Vec3d(x / radius_mid_surface, y / radius_mid_surface, z / radius_mid_surface);
				initializeSurfaceProperties(n_0, thickness);
			}
		}
	}
};
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody& body, const std::string& body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry() {};

private:
	void tagManually(size_t index_i)
	{
		SolidParticles& solid_particles = dynamic_cast<SolidParticles&>(base_particles_);
		if (solid_particles.pos0_[index_i][1] < 0.5 * particle_spacing_ref
			&& solid_particles.pos0_[index_i][1] > -0.5 * particle_spacing_ref)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/**
 *  The main program
 */
int main(int ac, char *av[])
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);

	/** Create a hemisphere body. */
	SolidBody hemisphere_body(system, makeShared<DefaultShape>("HemisphereBody"));
	hemisphere_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
	hemisphere_body.generateParticles<HemisphereParticleGenerator>();

	//hemisphere_body.addBodyStateForRecording<Vec3d>("PriorAcceleration");
	//hemisphere_body.addBodyStateForRecording<Vec3d>("NormalDirection");
	//hemisphere_body.addBodyStateForRecording<Vec3d>("InitialNormalDirection");
	//hemisphere_body.addBodyStateForRecording<Mat3d>("CorrectionMatrix");
	//hemisphere_body.addBodyStateForRecording<Mat3d>("TransformationMatrix");
	//hemisphere_body.addBodyStateForRecording<Mat3d>("DeformationGradient");
	//hemisphere_body.addBodyStateForRecording<Mat3d>("CorrectionMatrix");

	/** Define Observer. */
	ObserverBody hemisphere_observer(system, "HemisphereObserver");
	hemisphere_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerRelation hemisphere_body_inner(hemisphere_body);
	ContactRelation hemisphere_observer_contact(hemisphere_observer, { &hemisphere_body });

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	 /** Corrected configuration. */
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
		corrected_configuration(hemisphere_body_inner);
	/** Time step size calculation. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(hemisphere_body);
	/** stress relaxation. */
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
		stress_relaxation_first_half(hemisphere_body_inner, 3, true);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
		stress_relaxation_second_half(hemisphere_body_inner);
	SimpleDynamics<thin_structure_dynamics::DistributingPointForcesToShell>
		apply_point_forces(hemisphere_body, point_forces, reference_positions,
			time_to_full_external_force, particle_spacing_ref, 1.6);
	BoundaryGeometry boundary_geometry(hemisphere_body, "BoundaryGeometry");
	SimpleDynamics<solid_dynamics::FixedInAxisDirection> constrain_holder(boundary_geometry, Vecd(0.0, 1.0, 1.0));
	SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter>
		constrain_mass_center(hemisphere_body, Vecd(1.0, 1.0, 1.0));
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		hemisphere_position_damping(0.2, hemisphere_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		hemisphere_rotation_damping(0.2, hemisphere_body_inner, "AngularVelocity", physical_viscosity);
	/** Output */
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_hemisphere_max_displacement("Position", io_environment, hemisphere_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.exec();

	/**
	* From here the time stepping begines.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_hemisphere_max_displacement.writeToFile(0);
	
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 2.0;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < output_period)
		{
			if (ite % 100 == 0) {
				std::cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			apply_point_forces.exec(dt);

			stress_relaxation_first_half.exec(dt);
			constrain_mass_center.exec(dt);
			hemisphere_position_damping.exec(dt);
			hemisphere_rotation_damping.exec(dt);
			constrain_mass_center.exec(dt);
			stress_relaxation_second_half.exec(dt);

			ite++;
			dt = computing_time_step_size.exec();
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_hemisphere_max_displacement.writeToFile(ite);
		TickCount t2 = TickCount::now();
		write_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}

