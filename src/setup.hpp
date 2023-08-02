#pragma once

#include "defines.hpp"
#include "lbm.hpp"
#include "shapes.hpp"
#include "info.hpp"

#ifdef DEM
void main_setup(const LAMMPS_NS::LAMMPS* liggghts); // main setup script
#else
void main_setup(); // main setup script
#endif // DEM