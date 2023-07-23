/**
* @file 	dw_3d_pinched_semi-cylinder.cpp
* @brief 	This is the benchmark test of the shell.
* @details  We consider large deformation of a cylinder surface.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
* @version  0.1
 */
#include "sphinxsys.h"
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;
using namespace SPH;   // Namespace cite here.

/**
 * @brief Basic geometry parameters and numerical setup.
 */

class Parameter
{
protected:
	Real radius = 1.001;                                      /** Radius of the inner wall of the cylinder. */
	Real height = 3.048;                                     /** Height of the cylinder. */
	Real thickness = 0.03;                                    /** Thickness of the cylinder. */
	Real radius_mid_surface = radius + thickness / 2.0;      /** Radius of the mid surface. */
	int particle_number = 160;								 /** Particle number in the bottom edge. */
	Real particle_spacing_ref = 2.0 * Pi * radius_mid_surface / (Real)particle_number;
	int particle_number_height = height / particle_spacing_ref;
	int BWD = 1;
	Real BW = particle_spacing_ref * (Real)BWD;              /** Boundary width, determined by specific layer of boundary particles. */

	/** For material properties of the solid. */
	Real rho0_s = 1.0; 			                             /** Normalized density. */
	Real Youngs_modulus = 2.0685e7;	                         /** Normalized Youngs Modulus. */
	Real poisson = 0.3; 			                         /** Poisson ratio. */
	Real physical_viscosity = 2.0;                         /** physical damping, here we choose the same value as numerical viscosity. */

	std::vector<Vecd> reference_positions{ Vec3d(0.0, height, radius_mid_surface),
									   Vec3d(0.0, height, -radius_mid_surface) };
	Real time_to_full_external_force = 0.1;
};

/** Define application dependent particle generator for thin structure. */
class CylinderParticleGenerator : public SurfaceParticleGenerator, public Parameter
{
public:
	CylinderParticleGenerator(SPHBody& sph_body) : SurfaceParticleGenerator(sph_body) {};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number; i++)
		{
			for (int j = 0; j < (particle_number_height + BWD + 1); j++)
			{
				Real x = radius_mid_surface * cos(i * 2 * Pi / (Real)particle_number);
				Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
				Real z = radius_mid_surface * sin(i * 2 * Pi / (Real)particle_number);
				initializePositionAndVolumetricMeasure(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref);
				Vecd n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
				initializeSurfaceProperties(n_0, thickness);
			}
		}
	}
};
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle, public Parameter
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
		if (solid_particles.pos0_[index_i][1] < 0.0)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};
//----------------------------------------------------------------------
//  Define system, geometry, material, particles and all other things.
//----------------------------------------------------------------------
class PreSettingCase : public Parameter
{
protected:
	/** Domain bounds of the system. */
	BoundingBox system_domain_bounds;
	// Observer location
	StdVec<Vecd> observation_location = { Vec3d(0.0, height + 1.5 * particle_spacing_ref, radius_mid_surface) };
	/** Setup the system. */
	SPHSystem system;

	/** create a Cylinder body. */
	SolidBody cylinder_body;

	/** Define Observer. */
	ObserverBody cylinder_observer;

public:
	PreSettingCase() :
		system_domain_bounds(Vec3d(-radius - thickness, -BW, -radius - thickness),
			Vec3d(radius + thickness, height + BW, radius + thickness)),
		system(system_domain_bounds, particle_spacing_ref),
		cylinder_body(system, makeShared<DefaultShape>("CylinderBody")),
		cylinder_observer(system, "CylinderObserver")
	{
		//----------------------------------------------------------------------
		//	Creating bodies with corresponding materials and particles.
		//----------------------------------------------------------------------
		cylinder_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
		cylinder_body.generateParticles<CylinderParticleGenerator>();

		cylinder_observer.defineParticlesAndMaterial();
		cylinder_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	}
};
//----------------------------------------------------------------------
//  Define environment.
//----------------------------------------------------------------------
class Environment : public PreSettingCase
{
protected:
	Real loading_factor;
	/** Set body contact map
		 *  The contact map gives the data connections between the bodies
		 *  basically the the range of bodies to build neighbor particle lists
		 */
	InnerRelation cylinder_body_inner;
	ContactRelation cylinder_observer_contact;

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	 /** Corrected configuration. */
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
		corrected_configuration;
	/** Time step size calculation. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size;
	/** active-passive stress relaxation. */
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
		stress_relaxation_first_half;
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
		stress_relaxation_second_half; 
	BoundaryGeometry boundary_geometry; 
	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constrain_holder;
	SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter>
		constrain_mass_center;
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		cylinder_position_damping;
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		cylinder_rotation_damping;
	/** Output */
	IOEnvironment io_environment;
	BodyStatesRecordingToPlt write_states;
	ObservedQuantityRecording<Vecd> write_cylinder_max_displacement;
	/** Common particle dynamics. */
	std::vector<Vecd> point_forces;
	SimpleDynamics<thin_structure_dynamics::DistributingPointForcesToShell> apply_point_forces;


	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;

public:
	explicit Environment(Real loading_factor_new) :
		PreSettingCase(),
		loading_factor(loading_factor_new),
		cylinder_body_inner(cylinder_body),
		cylinder_observer_contact(cylinder_observer, { &cylinder_body }),
		corrected_configuration(cylinder_body_inner),
		computing_time_step_size(cylinder_body),
		stress_relaxation_first_half(cylinder_body_inner, 3, true),
		stress_relaxation_second_half(cylinder_body_inner), 
		boundary_geometry(cylinder_body, "BoundaryGeometry"),
		constrain_holder(boundary_geometry),
		constrain_mass_center(cylinder_body, Vecd(0.0, 0.0, 1.0)),
		cylinder_position_damping(0.2, cylinder_body_inner, "Velocity", physical_viscosity),
		cylinder_rotation_damping(0.2, cylinder_body_inner, "AngularVelocity", physical_viscosity),
		io_environment(system),
		write_states(io_environment, system.real_bodies_),
		write_cylinder_max_displacement("Position", io_environment, cylinder_observer_contact),
		point_forces({ Vec3d(0.0, 0.0, -loading_factor), Vec3d(0.0, 0.0, loading_factor) }),
		apply_point_forces(cylinder_body, point_forces, reference_positions,
			time_to_full_external_force, particle_spacing_ref, 1.6)
	{
		std::cout << "Running simulation for loading factor = " << loading_factor << "\n";
		/** Apply initial condition. */
		system.initializeSystemCellLinkedLists();
		system.initializeSystemConfigurations();
		corrected_configuration.exec();

		GlobalStaticVariables::physical_time_ = 0.0;
		write_states.writeToFile(0);
		write_cylinder_max_displacement.writeToFile(0);
	}

	virtual ~Environment() {};
	//----------------------------------------------------------------------
	//	For ctest.
	//----------------------------------------------------------------------
	int cmakeTest()
	{
		return 1;
	}

	/**
	 *  The main program
	 */
	void runCase()
	{

		/** Set the starting time.
		* From here the time stepping begins.
		*/
		GlobalStaticVariables::physical_time_ = 0.0;
		/** Setup physical parameters. */
		int ite = 0;
		Real end_time = 0.3;
		Real output_period = end_time / 100.0;
		Real dt = 0.0;

		/**
		 * Main loop
		 */
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
				apply_point_forces.exec(dt);

				stress_relaxation_first_half.exec(dt);
				constrain_holder.exec();
				constrain_mass_center.exec(dt);
				cylinder_position_damping.exec(dt);
				cylinder_rotation_damping.exec(dt);
				constrain_holder.exec();
				constrain_mass_center.exec(dt);
				stress_relaxation_second_half.exec(dt);

				ite++;
				dt = computing_time_step_size.exec();
				integral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_cylinder_max_displacement.writeToFile(ite);

			TickCount t2 = TickCount::now();
			//write_states.writeToFile();
			TickCount t3 = TickCount::now();
			interval += t3 - t2;
		}
		TickCount t4 = TickCount::now();

		write_states.writeToFile(1);

		TimeInterval tt;
		tt = t4 - t1 - interval;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	}
};

PYBIND11_MODULE(dw_3d_pinched_semi_cylinder_python, m)
{
	py::class_<Environment>(m, "pinched_semi_cylinder_from_sph_cpp")
		.def(py::init<const float&>())
		.def("CmakeTest", &Environment::cmakeTest)
		.def("RunCase", &Environment::runCase);
}

