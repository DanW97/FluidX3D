#pragma once

#include "defines.hpp"
#include "lbm.hpp"
#include "shapes.hpp"

#ifdef DEM
void main_setup(LAMMPS_NS::LAMMPS* lammps); // main setup script
#else
void main_setup(); // main setup script
#endif // DEM