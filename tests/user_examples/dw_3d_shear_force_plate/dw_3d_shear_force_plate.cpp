/**
* @file 	test_3d_thin_plate.cpp
* @brief 	This is the benchmark test of the shell.
* @details  We consider the body force applied on a quasi-static square plate.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
*/
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 10.0;                                          /** Length of the square plate. */
Real PW = 1.0;                                           /** Width of the square plate. */
Real PT = 0.1;                                           /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);                        /** Pseudo-normal. */
int particle_number = 5;								 /** Particle number in the direction of the width */
Real resolution_ref = PW / (Real)particle_number;        /** Initial reference particle spacing. */
int particle_number_PL = PL / resolution_ref;
int BWD = 3;
Real BW = resolution_ref * (Real)BWD;                    /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -0.5 * resolution_ref),
	Vec3d(PL + BW, PW + BW, 0.5 * resolution_ref));

/** For material properties of the solid. */
Real rho0_s = 1.0; 			                             /** Normalized density. */
Real Youngs_modulus = 1.2e6;	                         /** Normalized Youngs Modulus. */
Real poisson = 0.0; 			                         /** Poisson ratio. */
Real physical_viscosity = 20.0;                         /** physical damping, here we choose the same value as numerical viscosity. */

Real force = 1.2 * Youngs_modulus * PW * PT* PT* PT / 12.0 / PL / PL; /** Force applied on the plate end. */
Vec3d body_acceleration = Vec3d(0.0, 0.0, force / (particle_number * rho0_s * PT * resolution_ref * resolution_ref));
Real time_to_full_external_force = 1.0;

/** Define application dependent particle generator for thin structure. */
class PlateParticleGenerator : public ParticleGeneratorDirect
{
public:
	PlateParticleGenerator() : ParticleGeneratorDirect()
	{
		// the plate and boundary
		for (int i = 0; i < (particle_number_PL + BWD + 1); i++)
		{
			for (int j = 0; j < (particle_number); j++)
			{
				Real x = resolution_ref * i - BW;
				Real y = resolution_ref * j + 0.5 * resolution_ref;
				positions_volumes_.push_back(std::make_pair(Vecd(x, y, 0.0), resolution_ref * resolution_ref));
			}
		}
	}
};

/**
 * application dependent initial condition
 */
class PlateDynamicsInitialCondition
	: public thin_structure_dynamics::ShellDynamicsInitialCondition
{
public:
	PlateDynamicsInitialCondition(SolidBody &plate)
		: thin_structure_dynamics::ShellDynamicsInitialCondition(plate) {};
protected:
	void Update(size_t index_i, Real dt) override {
		/** initial pseudo-normal. */
		n_0_[index_i] = n_0;
		n_[index_i] = n_0;
		pseudo_n_[index_i] = n_0_[index_i];
	};
};

/** Define the geometry undergoing force. */
class GeometryUndergoingForce : public BodyPartByParticle
{
public:
	GeometryUndergoingForce(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&GeometryUndergoingForce::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~GeometryUndergoingForce() {};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_->pos_n_[index_i][0] > PL - 0.5 * resolution_ref)
		{
			body_part_particles_.push_back(index_i);
		}
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
	virtual ~BoundaryGeometry() {};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_->pos_n_[index_i][0] < 0.5 * resolution_ref)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public PartSimpleDynamicsByParticle
	, public DataDelegateSimple<ThinStructure, ShellParticles, ElasticSolid>
{

	StdLargeVec<Vecd> &dvel_dt_prior_;
public:
	TimeDependentExternalForce(SPHBody &sph_body, BodyPartByParticle &body_part)
		: PartSimpleDynamicsByParticle(sph_body, body_part),
		DataDelegateSimple<ThinStructure, ShellParticles, ElasticSolid>(sph_body),
		dvel_dt_prior_(particles_->dvel_dt_prior_) {}
	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		dvel_dt_prior_[index_i] = current_time < time_to_full_external_force ?
			current_time * body_acceleration / time_to_full_external_force : body_acceleration;
	}
};

/** Define application dependent observer particle generator. */
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(Vecd(PL + resolution_ref, 0.5 * PW, 0.0), 0.0));
	}
};

/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, resolution_ref);

	/** create a plate body. */
	ThinStructure plate_body(system, "PlateBody", makeShared<SPHAdaptation>(1.15, 1.0));
	/** create a plate material. */
	SharedPtr<LinearElasticSolid> plate_material = makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	/** create particles for the elastic body. */
	ShellParticles plate_body_particles(plate_body, plate_material, makeShared<PlateParticleGenerator>(), PT);
	plate_body_particles.addAVariableToWrite<Vec3d>("PriorAcceleration");

	/** Define Observer. */
	ObserverBody plate_observer(system, "PlateObserver");
	ObserverParticles observer_particles(plate_observer, makeShared<ObserverParticleGenerator>());

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	BodyRelationInner plate_body_inner(plate_body);
	BodyRelationContact plate_observer_contact(plate_observer, { &plate_body });

	/** Common particle dynamics. */
	GeometryUndergoingForce geometry_under_force(plate_body, "GeometryUndergoingForce");
	TimeDependentExternalForce impose_external_force(plate_body, geometry_under_force);

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	 /** initial condition */
	PlateDynamicsInitialCondition plate_initial_pseudo_normal(plate_body);
	 /** Corrected strong configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration_in_strong_form(plate_body_inner);
	/** Time step size calculation. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(plate_body);
	/** active-passive stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(plate_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(plate_body_inner);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
	thin_structure_dynamics::ConstrainShellBodyRegion constrain_holder(plate_body, boundary_geometry);
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>
		plate_position_damping(plate_body_inner, 0.2, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>
		plate_rotation_damping(plate_body_inner, 0.2, "AngularVelocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToPlt write_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<Vecd> write_plate_max_displacement("Position", in_output, plate_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	plate_initial_pseudo_normal.parallel_exec();
	corrected_configuration_in_strong_form.parallel_exec();

	/**
	* From here the time stepping begines.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_plate_max_displacement.writeToFile(0);
	
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 15.0;
	Real output_period = end_time / 200.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 100 == 0) {
				std::cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			impose_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			plate_position_damping.parallel_exec(dt);
			plate_rotation_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = 0.5 * computing_time_step_size.parallel_exec();
			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

		}
		write_plate_max_displacement.writeToFile(ite);
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}

