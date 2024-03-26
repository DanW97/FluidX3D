set -e
# command line argument(s): 
# -in [path to liggghts input file]
# -log [path to liggghts logs]
# device number (fluidx3d)

# TODO enforce that we at least get a liggghts file as a command arg

mkdir -p bin # create directory for executable
rm -f ./bin/LX3D # cleanup
rm -f ./post/*
# the only compiliation option here will be linux, no X11
bear -- mpiCC ./src/*.cpp -o ./bin/LX3D -std=c++17 -pthread -I./src/OpenCL/include -I$LIGGGHTS_PATH/include -L./src/OpenCL/lib -lOpenCL -L$LIGGGHTS_PATH/lib -l:libliggghts.so
# run LX3D
mpiexec -n 6 ./bin/LX3D
# ./bin/LX3D
