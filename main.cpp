/**
 * Sean McGlincy
 * Parallel Systems
 * Assignment 1
 *
 * Dependencies: I am using CLION as an IDE which uses CMAKE 3.8, and GCC, C+11
 * gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-11)
 * Running on Cento 7
 *
 *
 * Program:
 * The program begins by checking the terminal arguments before starting a high grade timer.
 * Five time measurements are taken and average together at the end of the program.
 * The program will read in a binary file of random intagers, positive or negative.
 *
 * The program copies these binary values to an array.  It is assumed that the file can be larger then the operating system
 * and the program will read in 1,000 integers at a time.  The program will also check to see if there are less then 1,000 integers
 * when importing the data.
 *
 * The first pass finds the min and max values.
 * An array is made with zero values.
 * The second pass increments the histogram.
 *  The program will print out to the terminal as well as create a text file after the timer has stopped.
 *
 * Running the File:
 * I recommend putting these in the same folder
 * To run file execute the binary file as ./filename
 * arg 1: is the filepath
 * arg 2: is the number of intervals
 * Example:  ./Assignment_01_Histogram random.binary 36
 *
 * */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <assert.h>
#include <mpi.h>

using namespace std;
using namespace chrono;

void clock(high_resolution_clock::time_point *array, int *time_samples) {
    for (int i = 0; i < *time_samples; i++) {
        array[i] = high_resolution_clock::now();
    }
}

double calculate_time(high_resolution_clock::time_point *start, high_resolution_clock::time_point *end, int *time_samples) {
    // Average time and convert to Micro Sec; 1 sec = 1,000,000 micro sec
    double total = 0;
    for (int i = 0; i < *time_samples; i++) {
        chrono::duration<double, std::milli> diff = end[i] - start[i];  // Time in Micro Sec
        total += diff.count();
    }
    return total / *time_samples;
}

//void print(int sample_size, int min, int max, int bucketSize, vector<int> *data) {
//    //////// Print //////////////
//    cout << "Sample Size: " << sample_size << endl;
//    cout << "Min Value: " << min << endl;
//    cout << "Max Value: " << max << endl;
//    cout << "bucket Range: " << bucketSize - 1 << endl;
//    for (unsigned int i = 0; i < data->size(); ++i) {
//        cout << "[" << min + (i * bucketSize) << ", " << min + ((i + 1) * bucketSize) - 1 << "] : " << &data[i] << endl;
//    }
//}
//
//int check_user_number(char *argv) {
//    char *endptr;
//    int intervalSize = strtol(argv, &endptr, 10);
//    if (!*argv || *endptr)
//        cerr << "Invalid number " << argv << '\n';
//    return intervalSize;
//}
//
//void broadcast() {}

void min_max_reduce(int *min, int *max, int *local_buffer, int &maxMessage, int root) {
        // todo make a function
}

void init_array(int* a, int arr_size){
    for(int i = 0; i< arr_size; i++)
        a[i] = 0;
}


int main(int argc, char *argv[]) {
//    cout << "Starting Program" << endl;


    //////// Start Clock //////////////
    // Use Chrono for high grade clock
    int time_samples = 5;
    high_resolution_clock::time_point clock_start[time_samples];
    high_resolution_clock::time_point clock_end[time_samples];
    clock(clock_start, &time_samples);

//    //////// USER INPUT//////////////
//    if (argc != 3) {
//        cout << "Error Error" << endl;
//        cout << "Please provide: binary data file and interval size" << endl;
//        exit(1);
//    }
//
//
//    assert(argc == 3);
//    string filePath = argv[1];
//    int intervalSize = check_user_number(argv[2]);
//    assert(intervalSize > 0);


//    int intervalSize = 10;
    string filePath = "random.binary";


    //////// MPI  Variables //////////////
    const int root = 0;
    int maxMessage = 2;
    int local_buffer[maxMessage];
    init_array(local_buffer, maxMessage);

    int local_min = numeric_limits<int>::max();
    int local_max = numeric_limits<int>::min();
    int num_iterations;
    int comm_sz;  // Number of process
    int my_rank;
    //////// MPI  INIT //////////////
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //////// Root Variables //////////////
    const int unit = sizeof(int);
    int bufferSize = unit * maxMessage * comm_sz;  // read in chuncks of file
    int size = bufferSize;
//    int bucketSize = 0;
    int fileLength = 0;
    int min = local_min;
    int max = local_max;
    int readBuffer[bufferSize];
    init_array(readBuffer, bufferSize);
    int seek = 0;





    if (my_rank == root) {
        cout << "Reading file"  << endl;
        //////// OPEN FILE //////////////
        ifstream fileInput;
        fileInput.open(filePath, ios::binary);
        if (fileInput.is_open()) {

            //////// Read FILE //////////////
            // Get file size
            fileInput.seekg(0, ios::end);
            fileLength = fileInput.tellg();
            fileInput.seekg(0);
            num_iterations = (int) ceil((double) fileLength / bufferSize);
            fileInput.close();
        } else {
            cout << "Can Not open file..." << endl;
            exit(1);
        }
    }

    MPI_Bcast(&num_iterations, 1, MPI_INT, root, MPI_COMM_WORLD );




    for (int i = 0; i < num_iterations; i++) {


        if (my_rank == 0) {
            cout << "reading File" << endl;
            //////// OPEN FILE //////////////
            ifstream fileInput;
            fileInput.open(filePath, ios::binary);
            if (fileInput.is_open()) {
                fileInput.seekg(seek);

                //////// First Pass:  Send to processors   //////////////
                //  Check if buffer is less then remainder of file
                if (fileLength - seek < size) {
                    size = fileLength - seek;   //todo update for mpi
                    maxMessage = size / comm_sz;
                }

                fileInput.read((char *) &readBuffer[0], size);  //todo fix


                seek = fileInput.tellg();
                fileInput.close();
            } else {
                cout << "Can Not open file..." << endl;
                exit(1);
            }
        }

        MPI_Bcast(&maxMessage, 1, MPI_INT, root, MPI_COMM_WORLD );
        MPI_Scatter(&readBuffer, maxMessage, MPI_INT, &local_buffer, maxMessage, MPI_INT, root, MPI_COMM_WORLD);

        for (int i = 0; i < maxMessage; i++) {
            if (local_min > local_buffer[i])
                local_min = local_buffer[i];
            if (local_max < local_buffer[i])
                local_max = local_buffer[i];
        }

        MPI_Reduce(&local_min, &min, 1, MPI_INT, MPI_MIN, root, MPI_COMM_WORLD);
        MPI_Reduce(&local_max, &max, 1, MPI_INT, MPI_MAX, root, MPI_COMM_WORLD);

        cout << i << " : Local Min: " << local_min  << " from processor " << my_rank << endl;
        cout << i << " : Local Max: " << local_max << " from processor " << my_rank  << endl;

    }


    ////////  END CLOCK //////////////
    //////// GET TIME //////////////

    if(my_rank == root){



    clock(clock_end, &time_samples);
    double total_time = calculate_time(clock_start, clock_end, &time_samples);
    cout << "AVG Time: " << total_time << " Milli Seconds" << endl;


        cout << "Program complete!" << endl;

        cout <<" : Min: " << min  << " from processor " << my_rank << endl;
        cout <<" : Max: " << max << " from processor " << my_rank  << endl;
    }

    MPI_Finalize();

    return 0;
}