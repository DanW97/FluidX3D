// LIGGGHTS
#include "lammps.h"
#include "input.h"
#include <mpi.h>

// FluidX3D
#include "setup.hpp"

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
    // liggghts -> input -> file(packing_file.c_str());
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
    sleep(60.0f);
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
    while(lbm.get_t() <= units.t(duration)) { // run for duration [s]
        // run liggghts for n steps
        const std::string cmd_string = "run "+std::to_string(couple_iter);
        const char *cmd = cmd_string.c_str();
        liggghts->input->one(cmd);
        // exchange forces
        if (lbm.get_t() % dump_iter == 0u) {
            // TODO add LIGGGHTS dumps too
            lbm.u.write_device_to_vtk(get_exe_path() + "../post/u_" + std::to_string(lbm.get_t()) + ".vtk");
            lbm.rho.write_device_to_vtk(get_exe_path() + "../post/rho_" + std::to_string(lbm.get_t()) + ".vtk");
        }
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
        impeller -> rotate(float3x3(float3(0.0f, 0.0f, 1.0f), domega)); // perform rotation
        lbm.voxelize_mesh_on_device(impeller, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));
        lbm.run(dt);
        // transfer to LIGGGHTS
        // wait for LIGGGHTS to run for n steps
    }
}


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