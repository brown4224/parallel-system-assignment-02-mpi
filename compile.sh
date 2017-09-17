#!bin/bash
rm ./mpi_program -f
mpiCC -g -Wall -o mpi_program main.cpp -std=c++0x 
mpiexec -n 5  mpi_program  ./random.binary  10



#scp '/home/sean/Desktop/School/parallel-systems/Assignment-02-MPI-Histogram/mpi_program' smcglin1@cloudland.kennesaw.edu:~/


