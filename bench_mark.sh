#! bin/bash
s1=2
s2=15


echo "Starting Bench Mark"
echo "Creating new Random number File"
/usr/bin/g++ -Wall  -o RandomNumber randomNumber.cpp -std=c++0x
./RandomNumber
echo "---------------------------------"
echo "Starting Base Line $s1"
./Assignment_01_Histogram ./random.binary $s1
echo "---------------------------------"
echo "Starting MPI Base Line $s1"
mpiexec -n 8  mpi_program  ./random.binary  $s1
echo "---------------------------------"
echo "Starting Base Line $s2"
./Assignment_01_Histogram ./random.binary $s2
echo "---------------------------------"
echo "Starting MPI Base Line $s2"
mpiexec -n 8  mpi_program  ./random.binary  $s2

