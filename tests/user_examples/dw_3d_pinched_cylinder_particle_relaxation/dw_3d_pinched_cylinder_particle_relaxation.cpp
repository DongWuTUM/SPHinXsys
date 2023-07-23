/**
 * @file 	test_3d_shell_particle_relaxation.cpp
 * @brief 	This is the test of using levelset to generate shell particles with single resolution and relax particles.
 * @details	We use this case to test the particle generation and relaxation by levelset for a complex thin structures geometry (3D).
 * @author 	Dong Wu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

Real radius_mid_surface = 5.0;                           /** Radius of the inner wall of the cylinder. */
Real height = 10.35;                                     /** Height of the cylinder. */
Real thickness_real = 0.094;                             /** Thickness of the cylinder. */

Real radius_large = radius_mid_surface + 0.5 * thickness_real;
Real radius_small = radius_mid_surface - 0.5 * thickness_real;
Real volumn_total = Pi * (radius_large * radius_large - radius_small * radius_small) * height;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_geometry = "./input/pullouted_cylinder.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 1.0e-3;
Vec3d domain_lower_bound(-5200.0 * scale, -5200.0 * scale, 0.0 * scale);
Vec3d domain_upper_bound(5200.0 * scale, 5200.0 * scale, 10350.0 * scale);
Vec3d cycle_center(-5100.0 * scale, -5100.0 * scale, 0.0 * scale);
Real dp_0 = 0.2;
Real thickness = 0.4;
// level set resolution much higher than that of particles is required
Real level_set_refinement_ratio = 5.0;
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

// Observer location
StdVec<Vecd> observation_location = { Vecd(radius_mid_surface, 0.0, 0.5 * height),
									  Vecd(0.0, radius_mid_surface, 0.5 * height),
									  Vecd(0.0, radius_mid_surface, 0.0) };

/** For material properties of the solid. */
Real rho0_s = 1.0; 			                             /** Normalized density. */
Real Youngs_modulus = 1.05e7;	                         /** Normalized Youngs Modulus. */
Real poisson = 0.3125; 			                         /** Poisson ratio. */
Real physical_viscosity = 2.0e2;                         /** physical damping, here we choose the same value as numerical viscosity. */

/** Difine point forces. */
Real F_full = 3.5e3;
std::vector<Vecd> point_forces{ Vec3d(F_full, 0.0, 0.0), Vec3d(-F_full, 0.0, 0.0) };
std::vector<Vecd> reference_positions{ Vec3d(radius_mid_surface, 0.0, height / 2.0),
									   Vec3d(-radius_mid_surface, 0.0, height / 2.0) };
Real time_to_full_external_force = 0.5;
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ImportedShellModel: public ComplexShape
{
public:
	explicit ImportedShellModel(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_geometry, cycle_center, scale);
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
		ShellParticles& shell_particles = dynamic_cast<ShellParticles&>(base_particles_);
		shell_particles.thickness_[index_i] = thickness_real;
	};
};
/**
 *  The main program
 */
int main(int ac, char* av[])
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, dp_0);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	system.setRunParticleRelaxation(true);
	/** Tag for starting with relaxed body-fitted particles distribution */
	system.setReloadParticles(true);
	#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
	#endif
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody cylinder_body(system, makeShared<ImportedShellModel>("ImportedShellModel"));
	cylinder_body.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	if (!system.RunParticleRelaxation() && system.ReloadParticles())
	{
		cylinder_body.generateParticles<ParticleGeneratorReload>(io_environment, cylinder_body.getName());
	}
	else
	{
		cylinder_body.defineBodyLevelSetShape(level_set_refinement_ratio);
		//here dummy linear elastic solid is use because no solid dynamics in particle relaxation
		cylinder_body.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness);
	}

	if (!system.RunParticleRelaxation() && !system.ReloadParticles())
	{
		std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
		return 0;
	}
	cylinder_body.addBodyStateForRecording<Vec3d>("PriorAcceleration");
	cylinder_body.addBodyStateForRecording<Vecd>("PseudoNormal");
	cylinder_body.addBodyStateForRecording<Vecd>("NormalDirection");
	cylinder_body.addBodyStateForRecording<Vec3d>("InitialNormalDirection");
	cylinder_body.addBodyStateForRecording<Real>("Thickness");

	/** Define Observer. */
	ObserverBody cylinder_observer(system, "CylinderObserver");
	cylinder_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerRelation cylinder_body_inner(cylinder_body);
	ContactRelation cylinder_observer_contact(cylinder_observer, { &cylinder_body });

	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	if (system.RunParticleRelaxation())
	{
		SimpleDynamics <RandomizeParticlePosition>  random_imported_model_particles(cylinder_body);
		// A  Physics relaxation step. 
		relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(cylinder_body_inner, thickness, level_set_refinement_ratio);
		relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(cylinder_body_inner, thickness, cos(0.1 * Pi));
		cylinder_body.addBodyStateForRecording<int>("UpdatedIndicator");
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_relaxed_particles(io_environment, system.real_bodies_);
		MeshRecordingToPlt write_mesh_cell_linked_list(io_environment, cylinder_body.getCellLinkedList());
		ReloadParticleIO write_particle_reload_files(io_environment, { &cylinder_body });
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_imported_model_particles.exec(0.25);
		relaxation_step_inner.mid_surface_bounding_.exec();
		write_relaxed_particles.writeToFile(0.0);
		cylinder_body.updateCellLinkedList();
		write_mesh_cell_linked_list.writeToFile(0.0);
		//----------------------------------------------------------------------
		//	Particle relaxation time stepping start here.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 2000)
		{
			if (ite_p % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_relaxed_particles.writeToFile(ite_p);
			}
			relaxation_step_inner.exec();
			ite_p += 1;
		}
		std::cout << "The physics relaxation process of imported model finish !" << std::endl;
		shell_normal_prediction.exec();
		write_relaxed_particles.writeToFile(ite_p);
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
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
	SimpleDynamics<thin_structure_dynamics::DistributingPointForcesToShell>
		apply_point_forces(cylinder_body, point_forces, reference_positions,
			time_to_full_external_force, dp_0, 2.3);

	BoundaryGeometry boundary_geometry(cylinder_body, "BoundaryGeometry");

	SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter>
		constrain_mass_center(cylinder_body, Vecd(1.0, 1.0, 1.0));
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		cylinder_position_damping(0.1, cylinder_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		cylinder_rotation_damping(0.1, cylinder_body_inner, "AngularVelocity", physical_viscosity);
	/** Output */
	BodyStatesRecordingToPlt write_states(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_cylinder_max_displacement("Position", io_environment, cylinder_observer_contact);

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
	write_cylinder_max_displacement.writeToFile(0);

	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 2.0;
	Real output_period = end_time / 200.0;
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
			cylinder_position_damping.exec(dt);
			cylinder_rotation_damping.exec(dt);
			constrain_mass_center.exec(dt);
			stress_relaxation_second_half.exec(dt);

			ite++;
			dt = 0.8 * computing_time_step_size.exec();
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_cylinder_max_displacement.writeToFile(ite);
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
