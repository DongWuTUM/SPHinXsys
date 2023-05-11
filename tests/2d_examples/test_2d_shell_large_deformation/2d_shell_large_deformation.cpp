/**
 * @file 	2d_shell.cpp
 * @brief 	This is the benchmark test of the 2d shell.
 * @details  We consider large deformation of a 2D shell.
 * @author 	Dong Wu, Chi Zhang and Xiangyu Hu
 * @version  0.1
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real radius = 24.5;											   /** Radius of the inner boundary of the cylinder. */
Real thickness = 1.0;										   /** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0;			   /** Radius of the mid surface. */
int particle_number = 2;									   /** Particle number in the thickness. */
Real particle_spacing_ref = thickness / (Real)particle_number; /** Initial reference particle spacing. */

// define arc length
Real deg = 215.0;
Real deg0 = -(deg-180)/2. / 180. * Pi;

int particle_number_mid_surface = 2 * radius_mid_surface * Pi * deg / 360.0 / particle_spacing_ref;

int BWD = 1;								/** number of boundary particles layers . */
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */

// for translating geometry from origin
Real move_x = 0.0;
Real move_y = 0.0;

BoundingBox system_domain_bounds(Vec2d(-radius - thickness + move_x, -radius - thickness + move_y),
								 Vec2d(radius + thickness + move_x, radius + thickness + move_y));
// Shell observer location
StdVec<Vecd> observation_location = {Vecd(0.0, radius_mid_surface)};
/** For material properties of the solid. */
Real rho0_s = 3.67346939;		  /** Normalized density. */
Real Youngs_modulus = 4.32e7;	  /** Normalized Youngs Modulus. */
Real poisson = 0.3;				  /** Poisson ratio. */
Real physical_viscosity = 2000.0; /** physical damping, here we choose the same value as numerical viscosity. */

Real time_to_full_external_force = 0.5;
Real gravitational_acceleration = -800.0;
/** Define application dependent particle generator for thin structure. */
class CylinderParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit CylinderParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		// the cylinder and boundary
		
		for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++)
		{
			Real x = radius_mid_surface *
					 cos(deg0 + (i + 0.5 - BWD) * deg / 360.0 * 2 * Pi / (Real)particle_number_mid_surface) + move_x;
			Real y = radius_mid_surface *
					 sin(deg0 + (i + 0.5 - BWD) * deg / 360.0 * 2 * Pi / (Real)particle_number_mid_surface) + move_y;
			initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_ref);
			Vec2d normal_direction = Vec2d((x-move_x) / radius_mid_surface, (y-move_y) / radius_mid_surface);
			initializeSurfaceProperties(normal_direction, thickness);
		}
	}
};

/** Define the displacement-control geometry. */
class DisControlGeometry : public BodyPartByParticle
{
public:
	DisControlGeometry(SPHBody& body, const std::string& body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&DisControlGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~DisControlGeometry() {};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_.pos_[index_i][0] < 0.5 * particle_spacing_ref && base_particles_.pos_[index_i][0] > -0.5 * particle_spacing_ref)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/** Define the controled displacement. */
class ControlDiaplacement : public thin_structure_dynamics::ConstrainShellBodyRegion
{
public:
	ControlDiaplacement(BodyPartByParticle& body_part)
		: ConstrainShellBodyRegion(body_part), 
		vel_(particles_->vel_){};
	virtual ~ControlDiaplacement() {};

protected:
	StdLargeVec<Vecd>& vel_;

	void update(size_t index_i, Real dt = 0.0)
	{
		vel_[index_i] = Vecd(0.0, -3.0);
	};
};

/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_.pos_[index_i][1] < radius_mid_surface * sin(deg0) + move_y)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity
{
public:
	explicit TimeDependentExternalForce(Vecd external_force)
		: Gravity(external_force) {}
	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_external_force ? current_time * global_acceleration_ / time_to_full_external_force : global_acceleration_;
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);

	/** Create a Cylinder body. */
	SolidBody cylinder_body(system, makeShared<DefaultShape>("CylinderBody"));
	cylinder_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
	cylinder_body.generateParticles<CylinderParticleGenerator>();

	// add more varaibles in output data
	cylinder_body.addBodyStateForRecording<Vecd>("PseudoNormal");
	cylinder_body.addBodyStateForRecording<Vecd>("Acceleration");
	cylinder_body.addBodyStateForRecording<Vecd>("AngularAcceleration"); 
	cylinder_body.addBodyStateForRecording<Vecd>("AngularVelocity");
	cylinder_body.addBodyStateForRecording<Matd>("DeformationGradient");
	cylinder_body.addBodyStateForRecording<Matd>("BendingDeformationGradient");
	cylinder_body.addBodyStateForRecording<Matd>("DeformationRate");
	cylinder_body.addBodyStateForRecording<Matd>("BendingDeformationGradientChangeRate");
	cylinder_body.addBodyStateForRecording<Matd>("GlobalStress");
	cylinder_body.addBodyStateForRecording<Matd>("GlobalMoment");
	cylinder_body.addBodyStateForRecording<Vecd>("GlobalShearStress");
	cylinder_body.addBodyStateForRecording<Vecd>("InitialNormalDirection");

	/** Define Observer. */
	ObserverBody cylinder_observer(system, "CylinderObserver");
	cylinder_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerRelation cylinder_body_inner(cylinder_body);
	ContactRelation cylinder_observer_contact(cylinder_observer, {&cylinder_body});

	/** Common particle dynamics. */
	SimpleDynamics<TimeStepInitialization> initialize_external_force(
		cylinder_body, makeShared<TimeDependentExternalForce>(Vec2d(0.0, gravitational_acceleration)));

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected configuration. */
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
		corrected_configuration(cylinder_body_inner);
	/** Time step size calculation. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(cylinder_body);
	/** stress relaxation. */
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
		stress_relaxation_first_half(cylinder_body_inner, 3, true);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
		stress_relaxation_second_half(cylinder_body_inner);
	/** Control the displacement. */
	DisControlGeometry dis_control_geometry(cylinder_body, "DisControlGeometry");
	SimpleDynamics<ControlDiaplacement> dis_control(dis_control_geometry);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(cylinder_body, "BoundaryGeometry");
	//SimpleDynamics<solid_dynamics::FixBodyPartConstraint> fixed_free_rotate_shell_boundary(boundary_geometry);
	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constrain_holder(boundary_geometry);
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
		cylinder_position_damping(0.2, cylinder_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
		cylinder_rotation_damping(0.2, cylinder_body_inner, "AngularVelocity", physical_viscosity);
	/** Output */
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_cylinder_max_displacement("Position", io_environment, cylinder_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.exec();

	/**
	 * From here the time stepping begins.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_cylinder_max_displacement.writeToFile(0);

	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 30.0;
	Real output_period = end_time / 1000.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	/**
	 * Main loop
	 */
	size_t numOutput = 0;
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			//initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.exec(dt);
			//fixed_free_rotate_shell_boundary.exec(dt);
			constrain_holder.exec(dt);
			dis_control.exec(dt);
			cylinder_position_damping.exec(dt);
			cylinder_rotation_damping.exec(dt);
			//fixed_free_rotate_shell_boundary.exec(dt);
			constrain_holder.exec(dt);
			dis_control.exec(dt);
			stress_relaxation_second_half.exec(dt);

			ite++;
			dt = computing_time_step_size.exec();

			//if (dt < 1.0e-8) {
			//	write_states.writeToFile();
			//}
			//if (ite > 236800) {
			//	write_states.writeToFile();
			//}

			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		//write_cylinder_max_displacement.writeToFile(ite);
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
