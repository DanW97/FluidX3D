// LIGGGHTS
#include "lammps.h"
#include "input.h"
#include <mpi.h>

// FluidX3D
#include "setup.hpp"

// Convenient dumping function
void dump_data(LBM &lbm, const std::string &path) {
    lbm.u.write_device_to_vtk(get_exe_path() + path + "/u_" + std::to_string(lbm.get_t()) + ".vtk");
    lbm.rho.write_device_to_vtk(get_exe_path() + path + "/rho_" + std::to_string(lbm.get_t()) + ".vtk");
    // lbm.dem_force.write_device_to_vtk(get_exe_path() + path + "/dem_force_" + std::to_string(lbm.get_t()) + ".vtk");
    // lbm.dem_torque.write_device_to_vtk(get_exe_path() + path + "/dem_torque_" + std::to_string(lbm.get_t()) + ".vtk");
}

/*void main_setup() { // Stokes drag validation; required extensions in defines.hpp: FORCE_FIELD, EQUILIBRIUM_BOUNDARIES
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint T = 100u; // check error every T steps
	const float R = 32.0f; // sphere radius
	const float Re = 0.01f; // Reynolds number
	const float nu = 1.0f; // kinematic shear viscosity
	const float rho = 1.0f; // density
	const uint L = to_uint(8.0f*R); // simulation box size
	const float u = units.u_from_Re(Re, 2.0f*R, nu); // velocity
	LBM lbm(L, L, L, nu); // flow driven by equilibrium boundaries
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E;
		if(sphere(x, y, z, lbm.center(), R)) {
			lbm.flags[n] = TYPE_S|TYPE_X; // flag boundary cells for force summation additionally with TYPE_X
		} else {
			lbm.rho[n] = units.rho_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R, rho, nu);
			const float3 un = units.u_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R);
			lbm.u.x[n] = un.x;
			lbm.u.y[n] = un.y;
			lbm.u.z[n] = un.z;
		}
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	double E1=1000.0, E2=1000.0;
	while(true) { // main simulation loop
		lbm.run(T);
		lbm.calculate_force_on_boundaries();
		lbm.F.read_from_device();
		const float3 force = lbm.calculate_force_on_object(TYPE_S|TYPE_X);
		const double F_theo = units.F_Stokes(rho, u, nu, R);
		const double F_sim = (double)length(force);
		const double E0 = fabs(F_sim-F_theo)/F_theo;
		print_info(to_string(lbm.get_t())+", expected: "+to_string(F_theo, 6u)+", measured: "+to_string(F_sim, 6u)+", error = "+to_string((float)(100.0*E0), 1u)+"%");
		if(converged(E2, E1, E0, 1E-4)) { // stop when error has sufficiently converged
			print_info("Error converged after "+to_string(lbm.get_t())+" steps to "+to_string(100.0*E0, 1u)+"%");
			wait();
			break;
		}
		E2 = E1;
		E1 = E0;
	}
} /**/

/*
void main_setup() { // mill - fluid only
    // dirty rounding lambda
    auto rnd = [](uint x) {return ((x + 5u) / 10u) * 10u;};
    // Simulation box size, viscosity, and other setups
    // mill is 321.24mm high, 152.18mm wide - aspect ratio = 2.11
    const uint Nx = 350u; const uint Ny = 350u; const uint Nz = 512u;
    const float eta = 0.1; // Pa s
    const float rho = 1.0f; // M
    const float si_rho = 1241.0; // kg m^-3
    const float si_nu = eta / si_rho; // m2 s^-1
    const float u = 0.2; // L T^-1
    const float si_Nx = 0.32124; // L
    const float si_u = 2.5; // m s^-1
    const uint liquid_height = Nz * 4u/5u;
    const uint stats_iter = 10u;
    units.set_m_kg_s(Nz, u, rho, si_Nx, si_u, si_rho);
    const float dump_frequency = 0.1; // every 100 ms
    const uint dump_iter = rnd(units.t(dump_frequency));
    // lbm object
    LBM lbm(Nx, Ny, Nz, units.nu(si_nu));
    // load in meshes
    const float chamber_size = units.x(0.32124f);
    const float impeller_size = units.x(0.29767f);
    const float3 center = lbm.center();
    // size parameter describes the *longest* side of the geometry bounding box (lbm units)
    Mesh* chamber = read_stl(get_exe_path() + "../mesh/mill.stl", lbm.size(), center, chamber_size);
    Mesh* impeller = read_stl(get_exe_path() + "../mesh/impeller.stl", lbm.size(), center, impeller_size);
    const float offset = 87.1f / 1000.0f; // height difference between impeller and chamber
    // impeller -> translate(float3(0.0f, 0.0f, units.x(offset)));
    const uint dt = 1; // every ~1ms
    const float si_omega = 2000.0 / 60.0 * 2.0 * 3.14; // 2000 rpm
    const float omega = units.omega(si_omega); // Lattice units
    const float domega = omega * dt;
    lbm.voxelize_mesh_on_device(chamber, TYPE_S, center);
    lbm.voxelize_mesh_on_device(impeller, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));

    // set initial conditions
    const ulong N=lbm.get_N(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
            if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
            if(lbm.flags[n]==0u) lbm.rho[n] = 1.;
            if(lbm.flags[n]==0u) lbm.u.x[n] = 0.;
            if(lbm.flags[n]==0u) lbm.u.y[n] = 0.;
            if(lbm.flags[n]==0u) lbm.u.z[n] = 0.;
    }

    // run sim
    lbm.run(0u);
    print_info("Lattice U  = "+to_string(u, 6u));
    print_info("Lattice Omega = "+to_string(units.omega(si_omega), 6u));
    print_info("Lattice Nu = "+to_string(units.nu_from_tau(lbm.get_tau()), 6u));
    print_info("Phys U  = "+to_string(si_u, 6u));
    print_info("Phys Omega = "+to_string(si_omega, 6u));
    print_info("Phys Nu = "+to_string(si_nu, 6u));
    print_info("Phys x (mm) = "+to_string(1000.0f * units.si_x(1.0f), 6u));
    print_info("Nz =  = "+to_string(units.x(chamber_size), 6u));
    print_info("1 LBM step = "+to_string(units.si_t(1u))+" s");
    print_info("File every "+to_string(rnd(units.t(dump_frequency)))+" LBM steps");
    print_info("dt "+to_string(dt)+" LBM steps");
    uint k = 0u;
    while(lbm.get_t() <= units.t(100.0f)) { // run for 100s
            // exchange forces
            if (lbm.get_t() % dump_iter == 0u) {
                    lbm.u.write_device_to_vtk(get_exe_path() + "../post/u_" + std::to_string(lbm.get_t()) + ".vtk");
                    lbm.rho.write_device_to_vtk(get_exe_path() + "../post/rho_" + std::to_string(lbm.get_t()) + ".vtk");
            }
            impeller -> rotate(float3x3(float3(0.0f, 0.0f, 1.0f), domega)); // perform rotation
            lbm.voxelize_mesh_on_device(impeller, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));
            lbm.run(dt);
    }
}
*/

/*
// Needs DEM define 
// TODO see which other defines are required
void main_setup(const LAMMPS_NS::LAMMPS* liggghts) { // mill - newtonian w/ solids
    // dirty rounding lambda
    auto rnd = [](uint x) {return ((x + 5u) / 10u) * 10u;};
    // Simulation box size, viscosity, and other setups
    // input file
    // TODO actually create this file
    const std::string packing_file = "../mill_packing.sim";
    // liggghts initialise
    liggghts -> input -> file(packing_file.c_str());
    const uint Nx = 350u; const uint Ny = 350u; const uint Nz = 512u;
    const float eta = 0.1; // Pa s
    const float rho = 1.0f; // M
    const float si_rho = 1241.0; // kg m^-3
    const float si_nu = eta / si_rho; // m2 s^-1
    const float u = 0.2; // L T^-1
    const float si_Nx = 0.32124; // L
    const float si_u = 2.5; // m s^-1
    const uint duration = 100; // s
    const uint liquid_height = Nz * 4u/5u;
    const uint stats_iter = 10u;
    units.set_m_kg_s(Nz, u, rho, si_Nx, si_u, si_rho);
    const float dump_frequency = 0.1; // every 100 ms
    const float couple_frequency = 0.01; // every 10 ms
    const uint dump_iter = rnd(units.t(dump_frequency));
    const uint couple_iter = rnd(units.t(couple_frequency));
    std::cout << "About to create LBM object" << std::endl;
    // lbm object
    LBM lbm(Nx, Ny, Nz, units.nu(si_nu), 0.0f, 0.0f, 0.0f, 5u, 20u);
    std::cout << "LBM object created" << std::endl;
    // important, gets all needed data from LIGGGHTS
    lbm.initialise_from_liggghts(liggghts);
    // load in meshes
    const float chamber_size = units.x(si_Nx);
    const float impeller_size = units.x(0.29767f);
    const float3 center = lbm.center();
    // size parameter describes the *longest* side of the geometry bounding box (lbm units)
    Mesh* chamber = read_stl(get_exe_path() + "../mesh/mill.stl", lbm.size(), center, chamber_size);
    Mesh* impeller = read_stl(get_exe_path() + "../mesh/impeller.stl", lbm.size(), center, impeller_size);
    const float offset = 87.1f / 1000.0f; // height difference between impeller and chamber
    // impeller -> translate(float3(0.0f, 0.0f, units.x(offset)));
    const uint dt = 1; // every ~1ms
    const float si_omega = 2000.0 / 60.0 * 2.0 * 3.14; // 2000 rpm
    const float omega = units.omega(si_omega); // Lattice units
    const float domega = omega * dt;
    std::cout << "Pre voxelise" << std::endl;
    lbm.voxelize_mesh_on_device(chamber, TYPE_S, center);
    lbm.voxelize_mesh_on_device(impeller, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));
    std::cout << "Post voxelise" << std::endl;
    // set initial conditions
    const ulong N=lbm.get_N(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
            if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
            if(lbm.flags[n]==0u) lbm.rho[n] = 1.;
            if(lbm.flags[n]==0u) lbm.u.x[n] = 0.;
            if(lbm.flags[n]==0u) lbm.u.y[n] = 0.;
            if(lbm.flags[n]==0u) lbm.u.z[n] = 0.;
    }

    // run sim
    // lbm initialise
    std::cout << "Pre GPU" << std::endl;
    lbm.run(0u);
    std::cout << "Post GPU" << std::endl;
    print_info("Lattice U  = "+to_string(u, 6u));
    print_info("Lattice Omega = "+to_string(units.omega(si_omega), 6u));
    print_info("Lattice Nu = "+to_string(units.nu_from_tau(lbm.get_tau()), 6u));
    print_info("Phys U  = "+to_string(si_u, 6u));
    print_info("Phys Omega = "+to_string(si_omega, 6u));
    print_info("Phys Nu = "+to_string(si_nu, 6u));
    print_info("Phys x (mm) = "+to_string(1000.0f * units.si_x(1.0f), 6u));
    print_info("Nz =  = "+to_string(units.x(chamber_size), 6u));
    print_info("1 LBM step = "+to_string(units.si_t(1u))+" s");
    print_info("File every "+to_string(rnd(units.t(dump_frequency)))+" LBM steps");
    print_info("dt "+to_string(dt)+" LBM steps");
    uint k = 0u;
    // inside this while loop, only execute fluidx3d on global rank 0
    // and include synchronisation barrier for coupling
    int rank; // rank in global communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    while(lbm.get_t() <= units.t(duration)) { // run for duration [s]
        // exchange forces
        if (rank == 0) {
            if (lbm.get_t() % couple_iter == 0u) {
                // get an idea of speed
                info.print_update();
                // reset hydrodynamic force and torque
                lbm.reset_coupling_forces();
                // send info to fluidx3d
                lbm.receive_liggghts_data(liggghts);
                // calculated lbm forces to particles
                lbm.send_force_data(liggghts);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // blocking call to ensure that force exchange has happened before continuing
        if (rank != 0) {
            // run liggghts for n steps
            const std::string cmd_string = "run "+std::to_string(couple_iter);
            const char *cmd = cmd_string.c_str();
            liggghts->input->one(cmd);
            // TODO add LIGGGHTS dumps too
        } else {
            // dump lbm data
            if (lbm.get_t() % dump_iter == 0u) {
                lbm.u.write_device_to_vtk(get_exe_path() + "../post/u_" + std::to_string(lbm.get_t()) + ".vtk");
                lbm.rho.write_device_to_vtk(get_exe_path() + "../post/rho_" + std::to_string(lbm.get_t()) + ".vtk");
            }
            impeller -> rotate(float3x3(float3(0.0f, 0.0f, 1.0f), domega)); // perform rotation
            lbm.voxelize_mesh_on_device(impeller, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));
            lbm.run(dt);
        }
    }
}
/**/

void main_setup(const LAMMPS_NS::LAMMPS* liggghts) { // flow over a sphere - validation of solid behaviour
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        auto rnd = [](uint x) {return ((x + 5u) / 10u) * 10u;};
        auto cd = [](float F, float u, float d, float rho) {return 2.0f * F / (rho * u * u * d * d * 0.5 * pi);};
        const float R = 16.0f; // sphere radius
        const float Re = 0.01f; // Reynolds number
        const float nu = 0.1f; // kinematic shear viscosity
        const float rho = 1.0f; // density
        const uint L = to_uint(16.0f*R); // simulation box size
        const float u = units.u_from_Re(Re, 2.0f*R, nu); // velocity
        // simulation control
        const uint dt = 1;
        const float dumping_frequency = 0.1; // 100ms
        const uint dump_iter = 10u;
        const uint stats_iter = 1u;

        // lbm object
        LBM lbm(L, L, L, nu, 0.0f, 0.0f, 0.0f, 1u, stats_iter); // flow driven by equilibrium boundaries
        lbm.initialise_from_liggghts(liggghts); // IMPORTANT TO CALL THIS

        // initial conditions
        const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
            if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
            if(!sphere(x, y, z, lbm.center(), R)) {
                lbm.rho[n] = units.rho_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R, rho, nu);
                const float3 un = units.u_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R);
                lbm.u.x[n] = un.x;
                lbm.u.y[n] = un.y;
                lbm.u.z[n] = un.z;
            }
        });
        // set sphere stuff
        float box_centre[3];
        box_centre[0] = lbm.center().x;
        box_centre[1] = lbm.center().y;
        box_centre[2] = lbm.center().z;
        for (uint i = 0u; i < 3u; i++) {
            lbm.dem_positions[i] = box_centre[i];
            // lbm.dem_positions[i] = 0.0f;
            lbm.dem_velocity[i] = 0.0f;
            lbm.dem_omega[i] = 0.0f;
        }
        lbm.dem_radii[0] = R;
        const float sqrt_r_half = sqrt(R * R - 0.5f);
        lbm.sphere_cap[0] = -R + 0.5f 
            + (1.0f / 12.0f - R * R)*atan(0.5f * sqrt_r_half / (0.5 - R * R))
            + (1.0f / 3.0f) * sqrt_r_half
            + (R * R - 1.0f / 12.0f)*atan(0.5f/sqrt_r_half)
            - (4.0f / 3.0f)*R*R*R*atan(0.25f / (sqrt_r_half * R));

        // run sim
        lbm.run(0u);
        print_info("Lattice U  = "+to_string(u, 6u));
        print_info("Lattice Nu = "+to_string(units.nu_from_tau(lbm.get_tau()), 6u));
        print_info("Phys U  = "+to_string(units.si_u(u), 6u));
        print_info("Phys Nu = "+to_string(units.si_nu(nu), 6u));
        print_info("Phys x (mm) = "+to_string(1000.0f * units.si_x(1.0f), 6u));
        print_info("Particle size (physical) = "+to_string(units.si_x(R*2.0f)));
        print_info("Particle size (lattice) = "+to_string(2.0f * R));
        print_info("Box centre (lattice) = "+to_string(lbm.center().x)+" "+to_string(lbm.center().y)+" "+to_string(lbm.center().z));
        print_info("1 LBM step = "+to_string(units.si_t(1u))+" s");
        print_info("File every "+to_string(rnd(units.t(dumping_frequency)))+" LBM steps");
        print_info("dt "+to_string(dt)+" LBM steps");
        while(lbm.get_t() <= 10u*stats_iter) { // run for 100s
            if (lbm.get_t() % stats_iter == 0) {
                lbm.dem_force.read_from_device();
                print_info("DEM force = "+to_string(-lbm.dem_force[0])+", "+to_string(lbm.dem_force[1])+", "+to_string(lbm.dem_force[2]));
                const float F = -lbm.dem_force[0];
                const float cd_theo = 24.0f / Re;
                const float cd_sim = cd(F, u, 2.0f*R, rho);
                const float percent_diff = 100.0f * fabs(cd_sim-cd_theo)/cd_theo;
                // print_info("t: "+to_string(lbm.get_t())+" cd_sim: "+to_string(cd_sim)+" cd_theo: "+to_string(cd_theo)+" percent difference: "+to_string(percent_diff)+"%");
                lbm.reset_coupling_forces();
            }
                
            // if (lbm.get_t() % dump_iter == 0u) dump_data(lbm, "../post");
            lbm.run(dt);
        }
    }
    // MPI_Barrier(MPI_COMM_WORLD);
} /**/

/*
// Needs DEM and a non-Newtonian define
void main_setup() { // mill - non-newtonian w/ solids
    // dirty rounding lambda
    auto rnd = [](uint x) {return ((x + 5u) / 10u) * 10u;};
    // Simulation box size, viscosity, and other setups
    // mill is 321.24mm high, 152.18mm wide - aspect ratio = 2.11
    const uint Nx = 350u; const uint Ny = 350u; const uint Nz = 512u;
    const float eta = 0.1; // Pa s
    const float rho = 1.0f; // M
    const float si_rho = 1241.0; // kg m^-3
    const float si_nu = eta / si_rho; // m2 s^-1
    const float u = 0.2; // L T^-1
    const float si_Nx = 0.32124; // L
    const float si_u = 2.5; // m s^-1
    const uint liquid_height = Nz * 4u/5u;
    const uint stats_iter = 10u;
    units.set_m_kg_s(Nz, u, rho, si_Nx, si_u, si_rho);
    const float dump_frequency = 0.1; // every 100 ms
    const uint dump_iter = rnd(units.t(dump_frequency));
    // lbm object
    LBM lbm(Nx, Ny, Nz, units.nu(si_nu));
    // load in meshes
    const float chamber_size = units.x(0.32124f);
    const float impeller_size = units.x(0.29767f);
    const float3 center = lbm.center();
    // size parameter describes the *longest* side of the geometry bounding box (lbm units)
    Mesh* chamber = read_stl(get_exe_path() + "../mesh/mill.stl", lbm.size(), center, chamber_size);
    Mesh* impeller = read_stl(get_exe_path() + "../mesh/impeller.stl", lbm.size(), center, impeller_size);
    const float offset = 87.1f / 1000.0f; // height difference between impeller and chamber
    // impeller -> translate(float3(0.0f, 0.0f, units.x(offset)));
    const uint dt = 1; // every ~1ms
    const float si_omega = 2000.0 / 60.0 * 2.0 * 3.14; // 2000 rpm
    const float omega = units.omega(si_omega); // Lattice units
    const float domega = omega * dt;
    lbm.voxelize_mesh_on_device(chamber, TYPE_S, center);
    lbm.voxelize_mesh_on_device(impeller, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));

    // set initial conditions
    const ulong N=lbm.get_N(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
        if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
        if(lbm.flags[n]==0u) lbm.rho[n] = 1.;
        if(lbm.flags[n]==0u) lbm.u.x[n] = 0.;
        if(lbm.flags[n]==0u) lbm.u.y[n] = 0.;
        if(lbm.flags[n]==0u) lbm.u.z[n] = 0.;
    }

    // run sim
    lbm.run(0u);
    print_info("Lattice U  = "+to_string(u, 6u));
    print_info("Lattice Omega = "+to_string(units.omega(si_omega), 6u));
    print_info("Lattice Nu = "+to_string(units.nu_from_tau(lbm.get_tau()), 6u));
    print_info("Phys U  = "+to_string(si_u, 6u));
    print_info("Phys Omega = "+to_string(si_omega, 6u));
    print_info("Phys Nu = "+to_string(si_nu, 6u));
    print_info("Phys x (mm) = "+to_string(1000.0f * units.si_x(1.0f), 6u));
    print_info("Nz =  = "+to_string(units.x(chamber_size), 6u));
    print_info("1 LBM step = "+to_string(units.si_t(1u))+" s");
    print_info("File every "+to_string(rnd(units.t(dump_frequency)))+" LBM steps");
    print_info("dt "+to_string(dt)+" LBM steps");
    uint k = 0u;
    while(lbm.get_t() <= units.t(100.0f)) { // run for 100s
        // exchange forces
        if (lbm.get_t() % dump_iter == 0u) {
            lbm.u.write_device_to_vtk(get_exe_path() + "../post/u_" + std::to_string(lbm.get_t()) + ".vtk");
            // lbm.rho.write_device_to_vtk(get_exe_path() + "../post/rho_" + std::to_string(lbm.get_t()) + ".vtk");
        }
        impeller -> rotate(float3x3(float3(0.0f, 0.0f, 1.0f), domega)); // perform rotation
        lbm.voxelize_mesh_on_device(impeller, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));
        lbm.run(dt);
        // transfer to LIGGGHTS
        // wait for LIGGGHTS to run for n steps
    }
}
*/