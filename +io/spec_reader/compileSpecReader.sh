module add gcc/7.3.0
module add hdf5_serial/1.10.5
module add openmpi/3.1.3

make init
make cleanAll
make
