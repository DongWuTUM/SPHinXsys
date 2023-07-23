/**
* @file 	dw_3d_pinched_hemisphere.cpp
* @brief 	This is the benchmark test of the shell.
* @details  We consider large deformation of a hemisphere surface.
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
	Real radius = 9.98;                                     /** Radius of the inner wall of the hemisphere. */
	Real thickness = 0.04;                                  /** Thickness of the hemisphere. */
	Real hole_angle = 18.0 / 180.0 * Pi;                    /** Angle of hole located at the hemisphere top. */
	Real radius_mid_surface = radius + thickness / 2.0;      /** Radius of the mid surface. */
	int particle_number = 80;								 /** Particle number in the bottom edge. */
	Real particle_spacing_ref = 2.0 * Pi * radius_mid_surface / (Real)particle_number;
	int particle_number_height = radius_mid_surface * (0.5 * Pi - hole_angle) / particle_spacing_ref;
	int BWD = 0;
	Real BW = particle_spacing_ref * (Real)BWD;              /** Boundary width, determined by specific layer of boundary particles. */

	/** For material properties of the solid. */
	Real rho0_s = 1.0; 			                             /** Normalized density. */
	Real Youngs_modulus = 6.825e7;	                         /** Normalized Youngs Modulus. */
	Real poisson = 0.3; 			                         /** Poisson ratio. */
	Real physical_viscosity = 60.0;                          /** physical damping, here we choose the same value as numerical viscosity. */

	std::vector<Vecd> reference_positions{ Vec3d(radius_mid_surface, 0.0, 0.0), Vec3d(-radius_mid_surface, 0.0, 0.0),
										   Vec3d(0.0, -radius_mid_surface, 0.0), Vec3d(0.0, radius_mid_surface, 0.0) };
	Real time_to_full_external_force = 0.1;
};

/** Define application dependent particle generator for thin structure. */
class HemisphereParticleGenerator : public SurfaceParticleGenerator, public Parameter
{
public:
	HemisphereParticleGenerator(SPHBody& sph_body) : SurfaceParticleGenerator(sph_body) {};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number_height + 1; i++)
		{
			Real circle_radius = radius_mid_surface * cos(i * (0.5 * Pi - hole_angle) / particle_number_height);
			int circle_particle_number = circle_radius * 2.0 * Pi / particle_spacing_ref;
			for (int j = 0; j < circle_particle_number; j++)
			{
				Real x = circle_radius * cos(j * 2 * Pi / (Real)circle_particle_number);
				Real y = circle_radius * sin(j * 2 * Pi / (Real)circle_particle_number);
				Real z = radius_mid_surface * sin(i * (0.5 * Pi - hole_angle) / particle_number_height);
				initializePositionAndVolumetricMeasure(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref);
				Vecd n_0 = Vec3d(x / radius_mid_surface, y / radius_mid_surface, z / radius_mid_surface);
				initializeSurfaceProperties(n_0, thickness);
			}
		}
	}
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
	StdVec<Vecd> observation_location = { Vecd(radius_mid_surface, 0.0, 0.0),
									  Vecd(0.0, radius_mid_surface, 0.0) };
	/** Setup the system. */
	SPHSystem system;

	/** create a hemisphere body. */
	SolidBody hemisphere_body;

	/** Define Observer. */
	ObserverBody hemisphere_observer; 

public:
	PreSettingCase() :
		system_domain_bounds(Vec3d(-radius - thickness, -radius - thickness, 0.0),
			Vec3d(radius + thickness, radius + thickness, radius + thickness)),
		system(system_domain_bounds, particle_spacing_ref),
		hemisphere_body(system, makeShared<DefaultShape>("HemisphereBody")),
		hemisphere_observer(system, "HemisphereObserver")
	{
		//----------------------------------------------------------------------
		//	Creating bodies with corresponding materials and particles.
		//----------------------------------------------------------------------
		hemisphere_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
		hemisphere_body.generateParticles<HemisphereParticleGenerator>();

		hemisphere_observer.defineParticlesAndMaterial();
		hemisphere_observer.generateParticles<ObserverParticleGenerator>(observation_location);
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
	InnerRelation hemisphere_body_inner;
	ContactRelation hemisphere_observer_contact;

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
		stress_relaxation_second_half; SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter>
		constrain_mass_center;
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		hemisphere_position_damping;
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		hemisphere_rotation_damping;
	/** Output */
	IOEnvironment io_environment;
	BodyStatesRecordingToPlt write_states;
	ObservedQuantityRecording<Vecd> write_hemisphere_max_displacement;
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
		hemisphere_body_inner(hemisphere_body),
		hemisphere_observer_contact(hemisphere_observer, { &hemisphere_body }),
		corrected_configuration(hemisphere_body_inner),
		computing_time_step_size(hemisphere_body),
		stress_relaxation_first_half(hemisphere_body_inner),
		stress_relaxation_second_half(hemisphere_body_inner),
		constrain_mass_center(hemisphere_body, Vecd(1.0, 1.0, 1.0)),
		hemisphere_position_damping(0.2, hemisphere_body_inner, "Velocity", physical_viscosity),
		hemisphere_rotation_damping(0.2, hemisphere_body_inner, "AngularVelocity", physical_viscosity),
		io_environment(system),
		write_states(io_environment, system.real_bodies_),
		write_hemisphere_max_displacement("Position", io_environment, hemisphere_observer_contact),
		point_forces({ Vec3d(loading_factor, 0.0, 0.0), Vec3d(-loading_factor, 0.0, 0.0),
										Vec3d(0.0, loading_factor, 0.0), Vec3d(0.0, -loading_factor, 0.0) }),
		apply_point_forces(hemisphere_body, point_forces, reference_positions,
			time_to_full_external_force, particle_spacing_ref, 1.15)
	{
		std::cout << "Running simulation for loading factor = " << loading_factor << "\n";
		/** Apply initial condition. */
		system.initializeSystemCellLinkedLists();
		system.initializeSystemConfigurations();
		corrected_configuration.exec();

		write_states.writeToFile(0);
		write_hemisphere_max_displacement.writeToFile(0);
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
		Real end_time = 2.0;
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
				constrain_mass_center.exec(dt);
				hemisphere_position_damping.exec(dt);
				hemisphere_rotation_damping.exec(dt);
				constrain_mass_center.exec(dt);
				stress_relaxation_second_half.exec(dt);

				ite++;
				dt = computing_time_step_size.exec();
				integral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_hemisphere_max_displacement.writeToFile(ite);

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

PYBIND11_MODULE(dw_3d_pinched_hemisphere_python, m)
{
	py::class_<Environment>(m, "pinched_hemisphere_from_sph_cpp")
		.def(py::init<const float&>())
		.def("CmakeTest", &Environment::cmakeTest)
		.def("RunCase", &Environment::runCase);
}

