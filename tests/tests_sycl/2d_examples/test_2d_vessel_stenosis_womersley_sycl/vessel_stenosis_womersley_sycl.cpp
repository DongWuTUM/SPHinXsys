/**
 * @file 	vessel_stenosis_womersley_sycl.cpp
 * @brief 	2D vessel stenosis with Womersley inflow example using SYCL.
 * @details This is a test case for Womersley flow in a 2D vessel with stenosis
 * @author 	Minhui Zhou, Dong Wu, Xiangyu Hu
 */
#include "sphinxsys.h" 
using namespace SPH;
//----------------------------------------------------------------------
// Select stenosis case (30% / 50% / 70%)
// By default we use the 30% stenosis CSVs.
// If you want 50% or 70% stenosis, change the filename suffix from "_0.3"
// to "_0.5" or "_0.7" accordingly (0.3=30%, 0.5=50%, 0.7=70%).
//----------------------------------------------------------------------
std::string womersley_velocity_profile_csv = "./input/womersley_velocity_profile_0.3.csv";
std::string outlet_pressure_csv = "./input/outlet_pressure_0.3.csv";
//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 0.0058;
Real DL1 = 4 * DH;
Real DL2 = 2 * DH;
Real DL3 = 16 * DH;
Real DL = DL1 + DL2 + DL3;
Real resolution_ref = DH / 12.0; 
Real BW = resolution_ref * 4.0; 
Real max_narrowing = 0.3; // chage to 0.5 or 0.7 for other stenosis cases
Real interpolationNum = 100;
BoundingBoxd system_domain_bounds(Vec2d(-DL1 - 0.5 * DL2 - BW, -0.5 * DH - BW), Vec2d(0.5 * DL2 + DL3, 0.5 * DH + BW));
//----------------------------------------------------------------------
//	Observation points
//----------------------------------------------------------------------
StdVec<Vecd> createAxialObservationPoints(Real full_length = DL, Vecd translation = Vecd(-DL1 - 0.5 * DL2, 0.0)) 
{
    StdVec<Vecd> observation_points;
    const int n_pts = 101;
    for (int i = 1; i < n_pts - 1; ++i) 
    {
        Real x = full_length * i / (n_pts - 1);
        Vecd point_coordinate(x, 0.0);
        observation_points.emplace_back(point_coordinate + translation);
    }
    return observation_points;
}
//----------------------------------------------------------------------
//	Material parameters of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1040.0;
Real mu_f = 0.004;
Real U_f = 4.44;
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Buffer parameters
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = Vec2d(-DL1 - 0.5 * DL2 + 2 * resolution_ref, 0);
Vec2d right_bidirectional_translation = Vec2d(0.5 * DL2 + DL3 - 2 * resolution_ref, 0);
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Inlet velocity profile for the left boundary
//----------------------------------------------------------------------
class WomersleyProfileCSV
{
  public:
    std::vector<Real> times;               
    std::vector<Real> radii;                  
    std::vector<std::vector<Real>> velocities; 
    Real period{0.0};

    explicit WomersleyProfileCSV(const std::string &csv_file, Real period_override = -1.0)
    {
        std::ifstream file(csv_file);
        if (!file.is_open())
        {
            std::cerr << "[ERROR] Cannot open file: " << csv_file << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::string line, token;
        std::getline(file, line);
        {
            std::stringstream ss(line);
            std::getline(ss, token, ','); 
            while (std::getline(ss, token, ','))
                radii.push_back(static_cast<Real>(std::stod(token)));
        }
        while (std::getline(file, line))
        {
            if (line.empty())
                continue;
            std::stringstream ls(line);
            std::getline(ls, token, ','); 
            if (token == "NaN" || token.empty())
                continue;
            times.push_back(static_cast<Real>(std::stod(token)));

            std::vector<Real> row;
            while (std::getline(ls, token, ','))
                row.push_back(static_cast<Real>(std::stod(token)));
            velocities.push_back(std::move(row));
        }
        if (period_override > 0.0)
        {
            period = period_override; 
        }
        else if (!times.empty())
        {
            period = times.back() - times.front();
        }

        std::cout << "[INFO] Loaded " << times.size()
                  << " time steps and " << radii.size()
                  << " radii. period = " << period << " s\n";
    }

    static size_t nearestIndex(const std::vector<Real> &v, Real x)
    {
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it == v.begin())
            return 0;
        if (it == v.end())
            return v.size() - 1;
        size_t hi = size_t(std::distance(v.begin(), it));
        size_t lo = hi - 1;
        return (std::fabs(x - v[lo]) <= std::fabs(x - v[hi])) ? lo : hi;
    }

    Real value(Real r, Real t) const
    {
        if (r <= radii.front())
            r = radii.front();
        else if (r >= radii.back())
            r = radii.back();
        if (period > 0.0 && !times.empty())
        {
            Real t0 = times.front();
            Real dt = t - t0;
            Real m = std::fmod(dt, period);
            if (m < Real(0))
                m += period;
            t = t0 + m;
        }

        size_t it = nearestIndex(times, t);
        size_t ir = nearestIndex(radii, r);
        return velocities[it][ir];
    }
};

class InflowVelocityPrescribed : public VelocityPrescribed<>
{
  public:
    WomersleyProfileCSV womersley_profile;
    
    InflowVelocityPrescribed(Real DH, Real U_f, Real mu_f,
                             const std::string &csv_file = womersley_velocity_profile_csv)
        : VelocityPrescribed<>(),
          womersley_profile(csv_file)
    {
        (void)DH;
        (void)U_f;
        (void)mu_f;
    }
    
    Real getAxisVelocity(const Vecd &input_position,
                         const Real &/*input_axis_velocity*/,
                         Real time)
    {
        Real r = std::fabs(input_position[1]);
        return womersley_profile.value(r, time);
    }
};

//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
class OutletPressureCSV
{
  public:
    std::vector<Real> times;
    std::vector<Real> pressures;
    Real period{0.0};

    explicit OutletPressureCSV(const std::string &csv_file, Real period_override = -1.0)
    {
        std::ifstream file(csv_file);
        if (!file.is_open())
        {
            std::cerr << "[ERROR] Cannot open file: " << csv_file << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::string line, tok_t, tok_p;

        while (std::getline(file, line))
        {
            if (line.empty())
                continue;
            std::stringstream ss(line);

            if (!std::getline(ss, tok_t, ','))
                continue;
            if (!std::getline(ss, tok_p, ','))
                continue;
            try
            {
                Real t = static_cast<Real>(std::stod(tok_t));
                Real p = static_cast<Real>(std::stod(tok_p));
                times.push_back(t);
                pressures.push_back(p);
            }
            catch (...)
            {
                continue;
            }
        }

        if (!times.empty())
        {
            if (period_override > Real(0))
            {
                period = period_override; 
            }
            else
            {
                period = times.back() - times.front(); 
            }
        }

        std::cout << "[INFO] Loaded " << times.size()
                  << " pressure samples. period = " << period << " s" << std::endl;
    }

    static size_t nearestIndex(const std::vector<Real> &v, Real x)
    {
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it == v.begin())
            return 0;
        if (it == v.end())
            return v.size() - 1;
        size_t hi = static_cast<size_t>(std::distance(v.begin(), it));
        size_t lo = hi - 1;
        return (std::fabs(x - v[lo]) <= std::fabs(x - v[hi])) ? lo : hi;
    }

    Real value(Real t) const
    {
        if (times.empty())
            return Real(0);

        if (period > Real(0))
        {
            Real t0 = times.front();
            Real dt = t - t0;
            Real m = std::fmod(dt, period);
            if (m < Real(0))
                m += period;
            t = t0 + m;
        }

        size_t it = nearestIndex(times, t);
        return pressures[it];
    }
};

class RightInflowPressure : public PressurePrescribed<>
{
public:
    OutletPressureCSV pressure_profile;

    RightInflowPressure(Real /*a para*/,
                        const std::string &csv_file = outlet_pressure_csv)
        : PressurePrescribed<>(0.0),
          pressure_profile(csv_file)
    {}

    Real getPressure(const Real& /*input_pressure*/, Real time)
    {
        return 1000.0 * pressure_profile.value(time); 
    }
};
//----------------------------------------------------------------------
//	stenosis definition
//  a - maximaum narrowing, X0 - halflength of stenosis
//----------------------------------------------------------------------

Real outline(Real x_rel, Real a, Real X0)
{
    if (std::abs(x_rel) <= X0)
    {
        return 1.0 - 0.5 * a * (1.0 + std::cos(Pi * x_rel / X0));
    }
    else
    {
        return 1.0;
    }
}

std::vector<Vecd> createStenosisUpper(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    Real dx = (2.0 * X0) / N; 
    std::vector<Vecd> stenosis_upper;
    stenosis_upper.reserve(N + 1);

    for (int i = N; i >= 0; --i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = (D * 0.5) * y_rel;
        stenosis_upper.push_back(Vecd(x_rel, +y));
    }
    stenosis_upper.push_back(stenosis_upper.front());
    return stenosis_upper;
}

std::vector<Vecd> createStenosisLower(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    Real dx = (2.0 * X0) / N;
    std::vector<Vecd> stenosis_lower;
    stenosis_lower.reserve(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = -(D * 0.5) * y_rel;
        stenosis_lower.push_back(Vecd(x_rel, y));
    }
    stenosis_lower.push_back(stenosis_lower.front());
    return stenosis_lower;
}
//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
std::vector<Vecd> createBloodShape()
{
    // geometry
    std::vector<Vecd> blood_shape;
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH));
    blood_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH));
    blood_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH));
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    return blood_shape;
}
class Blood : public ComplexShape
{
  public:
    explicit Blood(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon blood_outer_boundary(createBloodShape());
        add<MultiPolygonShape>(blood_outer_boundary, "BloodOuterBoundary");

        MultiPolygon stenosisUpper(createStenosisUpper(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(stenosisUpper);
        MultiPolygon stenosisLower(createStenosisLower(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(stenosisLower);
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH + BW));
    outer_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH + BW));
    outer_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH - BW));
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));
    return outer_wall_shape;
}
std::vector<Vecd> createCompleteInnerWallShape(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    std::vector<Vecd> lumen;

    lumen.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * D));
    lumen.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * D));
    lumen.push_back(Vecd(-X0, 0.5 * D));

    Real dx = (2.0 * X0) / N;
    std::vector<Vecd> stenosis_lower1;
    stenosis_lower1.reserve(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = (D * 0.5) * y_rel;
        lumen.push_back(Vecd(x_rel, +y));
    }

    lumen.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * D));

    lumen.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * D));
    lumen.push_back(Vecd(X0, -0.5 * D));

    for (int i = N; i >= 0; --i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = -(D * 0.5) * y_rel;
        lumen.push_back(Vecd(x_rel, y));
    }

    lumen.push_back(lumen.front());
    return lumen;
}
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon boundary_outer_wall_shape(createOuterWallShape());
        add<MultiPolygonShape>(boundary_outer_wall_shape, "OuterWallBoundary");
        MultiPolygon boundary_inner_wall_shape(createCompleteInnerWallShape(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(boundary_inner_wall_shape, "OuterWallBoundary");
    }
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // sph_system.setRunParticleRelaxation(true);
    // sph_system.setReloadParticles(false);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.setGenerateRegressionData(true);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody blood(sph_system, makeShared<Blood>("WaterBody"));
    blood.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> particle_buffer(4.0);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();

    ObserverBody velocity_axial_observer(sph_system, "VelocityAxialObserver");
    velocity_axial_observer.defineAdaptationRatios(0.25, 1.0); 
    velocity_axial_observer.generateParticles<ObserverParticles>(createAxialObservationPoints());
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {  
        LevelSetShape *blood_level_set_shape = blood.defineBodyLevelSetShape(2.0)
                                                ->correctLevelSetSign()
                                                ->cleanLevelSet()
                                                ->writeLevelSet();
        blood.generateParticlesWithReserve<BaseParticles, Lattice>(particle_buffer);

        LevelSetShape *wall_boundary_level_set_shape = wall_boundary.defineBodyLevelSetShape(2.0)
                                                ->correctLevelSetSign()
                                                ->cleanLevelSet()
                                                ->writeLevelSet();
        wall_boundary.generateParticles<BaseParticles, Lattice>();
        //----------------------------------------------------------------------
        //	Creating body parts.
        //----------------------------------------------------------------------
        NearShapeSurface near_blood_surface(blood);
        NearShapeSurface near_wall_boundary_surface(wall_boundary);
        //----------------------------------------------------------------------
        //  Define body relation map.
        //----------------------------------------------------------------------
        Inner<> blood_inner(blood);
        Inner<> wall_boundary_inner(wall_boundary);
        Contact<> blood_contact(blood, {&wall_boundary});
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        SPHSolver sph_solver(sph_system);
        auto &main_methods = sph_solver.addParticleMethodContainer(par_host);
        //----------------------------------------------------------------------
        // Define the numerical methods used in the simulation.
        // Note that there may be data dependence on the sequence of constructions.
        // Generally, the configuration dynamics, such as update cell linked list,
        // update body relations, are defiend first.
        // Then the geometric models or simple objects without data dependencies,
        // such as gravity, initialized normal direction.
        // After that, the major physical particle dynamics model should be introduced.
        // Finally, the auxiliary models such as time step estimator, initial condition,
        // boundary condition and other constraints should be defined.
        //----------------------------------------------------------------------
        auto &wall_boundary_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
        auto &blood_cell_linked_list = main_methods.addCellLinkedListDynamics(blood);
        auto &wall_boundary_update_inner_relation = main_methods.addRelationDynamics(wall_boundary_inner);
        auto &blood_update_complex_relation = main_methods.addRelationDynamics(blood_inner, blood_contact);

        auto &random_wall_boundary_particles = main_methods.addStateDynamics<RandomizeParticlePositionCK>(wall_boundary);
        auto &random_blood_particles = main_methods.addStateDynamics<RandomizeParticlePositionCK>(blood);

        auto &wall_boundary_relaxation_residual =
            main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(wall_boundary_inner)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(wall_boundary, *wall_boundary_level_set_shape);
        auto &wall_boundary_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(wall_boundary);
        auto &wall_boundary_update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(wall_boundary);
        auto &wall_boundary_level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_wall_boundary_surface);

        auto &blood_relaxation_residual =
            main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(blood_inner)
                .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(blood_contact)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(blood, *blood_level_set_shape);
        auto &blood_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(blood);
        auto &blood_update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(blood);
        auto &blood_level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_blood_surface);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
        auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&wall_boundary, &blood});
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        random_wall_boundary_particles.exec();
        random_blood_particles.exec();
        //----------------------------------------------------------------------
        //	First output before the simulation.
        //----------------------------------------------------------------------
        body_state_recorder.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            wall_boundary_cell_linked_list.exec();
            wall_boundary_update_inner_relation.exec();

            wall_boundary_relaxation_residual.exec();
            Real wall_boundary_relaxation_step = wall_boundary_relaxation_scaling.exec();
            wall_boundary_update_particle_position.exec(wall_boundary_relaxation_step);
            wall_boundary_level_set_bounding.exec();

            blood_cell_linked_list.exec();
            blood_update_complex_relation.exec();
            blood_relaxation_residual.exec();
            Real blood_relaxation_step = blood_relaxation_scaling.exec();
            blood_update_particle_position.exec(blood_relaxation_step);
            blood_level_set_bounding.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;
        write_particle_reload_files.writeToFile();
        return 0;
    }

    //----------------------------------------------------------------------
    //	Particle reloading.
    //----------------------------------------------------------------------
    blood.generateParticlesWithReserve<BaseParticles, Reload>(particle_buffer, blood.getName());
    wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName());
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    AlignedBox left_emitter_shape(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize);
    AlignedBoxByCell left_emitter(blood, left_emitter_shape);

    AlignedBox right_emitter_shape(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize);
    AlignedBoxByCell right_emitter(blood, right_emitter_shape);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    // ----------------------------------------------------------------------
    Inner<> blood_inner(blood);
    Inner<> wall_boundary_inner(wall_boundary);
    Contact<> blood_wall_contact(blood, {&wall_boundary});
    Contact<> velocity_observer_contact_axial(velocity_axial_observer, {&blood});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> blood_cell_linked_list(blood);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> blood_body_update_complex_relation(blood_inner, blood_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> fluid_observer_contact_relation(velocity_observer_contact_axial);
    ParticleSortCK<MainExecutionPolicy> particle_sort(blood);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal_direction(wall_boundary); // run on CPU
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> blood_advection_step_setup(blood);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> blood_update_particle_position(blood);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(blood_inner, 0.5), blood_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(blood_inner, blood_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallNoRiemannCK>
        fluid_acoustic_step_2nd_half(blood_inner, blood_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplexInternalPressureBoundary>
        fluid_density_regularization(blood_inner, blood_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(blood_inner, blood_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionComplexBulkParticlesCK>
        transport_correction_ck(blood_inner, blood_wall_contact);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(blood, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> fluid_acoustic_time_step(blood);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(blood_inner, blood_wall_contact);
    fluid_dynamics::BidirectionalBoundaryCK<MainExecutionPolicy, LinearCorrectionCK, InflowVelocityPrescribed>
        bidirectional_velocity_condition_left(left_emitter, DH, U_f, mu_f);
    fluid_dynamics::BidirectionalBoundaryCK<MainExecutionPolicy, LinearCorrectionCK, RightInflowPressure>
        bidirectional_pressure_condition_right(right_emitter, 0.0);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::OutflowParticleDeletion> out_flow_particle_deletion(blood);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(blood, "Pressure");
    body_states_recording.addToWrite<int>(blood, "BufferIndicator");
    RestartIOCK<MainExecutionPolicy> restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<MainExecutionPolicy, Vecd, RestoringCorrection>> write_centerline_velocity("Velocity", velocity_observer_contact_axial);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        sv_physical_time->setValue(restart_io.readRestartFiles(sph_system.RestartStep()));
    }
    wall_normal_direction.exec();
    blood_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    blood_body_update_complex_relation.exec();
    fluid_observer_contact_relation.exec();
    fluid_boundary_indicator.exec();
    bidirectional_velocity_condition_left.tagBufferParticles();
    bidirectional_pressure_condition_right.tagBufferParticles();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    size_t screen_output_interval = 100;
    size_t observation_sample_interval = screen_output_interval * 2;
    size_t restart_output_interval = screen_output_interval * 10;
    Real end_time = 1.6; // T = 0.8s
    Real output_interval = 0.01;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_centerline_velocity.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            time_instance = TickCount::now();
            fluid_density_regularization.exec();
            blood_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            transport_correction_ck.exec();
            fluid_viscous_force.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                acoustic_dt = SMIN(fluid_acoustic_time_step.exec(), advection_dt);
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                bidirectional_velocity_condition_left.applyBoundaryCondition(acoustic_dt);
                bidirectional_pressure_condition_right.applyBoundaryCondition(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            blood_update_particle_position.exec();
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body observables and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerline_velocity.writeToFile(number_of_iterations);
                }

                if (number_of_iterations % restart_output_interval == 0)
                {
                    restart_io.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            time_instance = TickCount::now();
            
            /** inflow emitter injection*/
            bidirectional_velocity_condition_left.injectParticles();
            bidirectional_pressure_condition_right.injectParticles();
            bidirectional_velocity_condition_left.indicateOutFlowParticles();
            bidirectional_pressure_condition_right.indicateOutFlowParticles();
            out_flow_particle_deletion.exec();
            /** Update cell linked list and configuration. */
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sort.exec();
            }
            blood_cell_linked_list.exec();
            blood_body_update_complex_relation.exec();
            fluid_observer_contact_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
            fluid_boundary_indicator.exec();
            bidirectional_velocity_condition_left.tagBufferParticles();
            bidirectional_pressure_condition_right.tagBufferParticles();
        }

        TickCount t2 = TickCount::now();

        body_states_recording.writeToFile();
        fluid_observer_contact_relation.exec();

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
    //----------------------------------------------------------------------
    // Post-run regression test to ensure that the case is validated
    // ----------------------------------------------------------------------
    if (sph_system.GenerateRegressionData())
    {
        write_centerline_velocity.generateDataBase(0.05);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_centerline_velocity.testResult();
    }

    return 0;
}
