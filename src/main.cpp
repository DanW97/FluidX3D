#pragma once
// LIGGGHTS
#include "lammps.h"
#include <mpi.h>
#include "input.h"
#include <string.h>
#include <signal.h>
#include "signal_handling.h"

// FluidX3D
#include "info.hpp"
#include "lbm.hpp"
#include "setup.hpp"

#ifdef DEM
void main_physics(LAMMPS_NS::LAMMPS* lammps) {
        info.print_logo();
        main_setup(lammps); // execute setup
        running = false;
        exit(0); // make sure that the program stops
}
int main(int argc, char *argv[]) {
    // TODO separate args for liggghts and fluidx3d (device id)
    struct sigaction int_action, usr1_action, term_action;
    memset(&int_action, 0, sizeof(struct sigaction));
    memset(&usr1_action, 0, sizeof(struct sigaction));
    memset(&term_action, 0, sizeof(struct sigaction));

    int_action.sa_handler = LAMMPS_NS::SignalHandler::int_handler;
    sigaction(SIGINT, &int_action, NULL);
    // SIGTERM is handled the same way as sigint. Note that OpenMPI (and possibly other flavours)
    // convert a SIGINT to mpirun to a SIGTERM to its children. That's why we need to catch it too.
    sigaction(SIGTERM, &int_action, NULL);
    usr1_action.sa_handler = LAMMPS_NS::SignalHandler::usr1_handler;
    sigaction(SIGUSR1, &usr1_action, NULL);

    // read in LIGGGHTS file and create liggghts object
    MPI_Init(&argc, &argv);
    const LAMMPS_NS::LAMMPS *lammps = new LAMMPS_NS::LAMMPS(argc, argv, MPI_COMM_WORLD);
    // main_arguments = get_main_arguments(argc, argv);
    thread compute_thread(main_physics, lammps);
    do { // main console loop
            info.print_update();
            sleep(0.050);
    } while(running);
    compute_thread.join();
    return 0;
}
#else
void main_physics() {
        info.print_logo();
        main_setup(); // execute setup
        running = false;
        exit(0); // make sure that the program stops
}
int main(int argc, char* argv[]) {
        main_arguments = get_main_arguments(argc, argv);
        thread compute_thread(main_physics());
        do { // main console loop
                info.print_update();
                sleep(0.050);
        } while(running);
        compute_thread.join();
        return 0;
}
#endif // DEM

