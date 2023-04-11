set -e
# command line argument(s): 
# -in [path to liggghts input file]
# -log [path to liggghts logs]
# device number (fluidx3d)

# TODO enforce that we at least get a liggghts file as a command arg

mkdir -p bin # create directory for executable
rm -f ./bin/LX3D # cleanup

# the only compiliation option here will be linux, no X11
mpiCC ./src/*.cpp -o ./bin/LX3D -std=c++17 -pthread -I./src/OpenCL/include -I$LIGGGHTS_DIR/include -L./src/OpenCL/lib -lOpenCL -L/home/dan/LIGGGHTS/build -l:libliggghts.so

./bin/LX3D "$@" # run LX3D
