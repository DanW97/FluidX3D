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

#include "coupling.hpp"

// do setup in here
void lbm_setup() {
	info.print_logo();
	// TODO 
	// running = false;
	exit(0); // make sure that the program stops
}

// write in here what constitutes an lbm iteration
void lbm_iteration(int n_iter) {

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

    LAMMPS_NS::LAMMPS *lammps = new LAMMPS_NS::LAMMPS(argc, argv, MPI_COMM_WORLD);
    // LBM setup and create lbm object

    // TODO ensure that only 1 host rank uses fluidx3d
    // initial iteration
    thread fx3d_thread(lbm_setup);
	info.print_update();
    lammps -> input -> file();
    // blocking call for fx3d
	fx3d_thread.join();

    //  particle positions to coverage collision operator

    //  calculated forces to particles

    // TODO have some nice definition for this that allows
    // lbm and dem to sync
    const uint n_iters = 1000;
    const uint n_iter = 10;
    // main loop
    uint iter = 0;
    while (iter <= n_iters) {
        // lbm
        // TODO see if the thread needs to be wrapped up in mpi or not
        thread fx3d_thread(lbm_iteration, n_iter);
        // dem
        const std::string cmd_string = "run "+std::to_string(n_iter);
        const char *cmd = cmd_string.c_str();
        lammps->input->one(cmd);

        // particle positions to coverage collision operator

        // calculated lbm forces to particles

        // data output?
        // TODO ensure that lbm forces can be dumped separately to the contact ones


        iter += n_iter;
    }


	return 0;
}