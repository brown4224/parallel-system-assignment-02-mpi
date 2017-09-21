/**
 * Sean McGlincy
 * Parallel Systems
 * Assignment 2
 *
 * Dependencies: I am using CLION as an IDE which uses CMAKE 3.8, and GCC, C+11
 * gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-11)
 * Running on Cento 7
 *
 *
 * Program:
 * Running the File:
 * I recommend putting these in the same folder
 * To run file execute the binary file as ./filename
 * arg 1: is the filepath
 * arg 2: is the number of intervals
 * Example: mpiCC -g -Wall -o mpi_program main.cpp -std=c++0x
 *          mpiexec -n 5  mpi_program  ./random.binary  10
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

double
calculate_time(high_resolution_clock::time_point *start, high_resolution_clock::time_point *end, int *time_samples) {
    // Average time and convert to Micro Sec; 1 sec = 1,000,000 micro sec
    double total = 0;
    for (int i = 0; i < *time_samples; i++) {
        chrono::duration<double, std::milli> diff = end[i] - start[i];  // Time in Micro Sec
        total += diff.count();
    }
    return total / *time_samples;
}

void print(int sample_size, int min, int max, int bucketSize, int *data, int myrank, int intervalSize) {
    cout << "Report from processor:  " << myrank << endl;
    cout << "Sample Size: " << sample_size << endl;
    cout << "Min Value: " << min << endl;
    cout << "Max Value: " << max << endl;
    cout << "bucket Range: " << bucketSize - 1 << endl;
    for (int i = 0; i < intervalSize; i++) {
        cout << "[" << min + (i * bucketSize) << ", " << min + ((i + 1) * bucketSize) - 1 << "] : " << data[i] << endl;
    }
}

int check_user_number(char *argv) {
    char *endptr;
    int intervalSize = strtol(argv, &endptr, 10);
    if (!*argv || *endptr)
        cerr << "Invalid number " << argv << '\n';
    return intervalSize;
}


void init_array(int *a, int arr_size) {
    new int[arr_size];
    for (int i = 0; i < arr_size; i++)
        a[i] = 0;
}


void readFile(string *filePath, int *readBuffer, int *messageSize, int *fileLength, int *comm_sz, int *seek, int *size,
              int *unit, int *remainder) {
    //////// OPEN FILE //////////////
    ifstream fileInput;
    fileInput.open(*filePath, ios::binary);
    if (fileInput.is_open()) {
        fileInput.seekg(*seek);

        //  Check if buffer is less then remainder of file
        if (*fileLength - *seek < *size) {
            *size = *fileLength - *seek;
            *messageSize = (*size / *unit) / *comm_sz;
            *remainder = (*size / *unit) - (*messageSize * *comm_sz);
        }
        fileInput.read((char *) &readBuffer[0], *size);
        *seek = fileInput.tellg();
        fileInput.close();
    } else {
        cout << "Can Not open file..." << endl;
        exit(1);
    }
}


void send_data(int *readBuffer, int *local_buffer, int *messageSize, int root, int *my_rank, int *size, int *unit,
               int *remainder) {
//    MPI_Bcast(messageSize, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Scatter(readBuffer, *messageSize, MPI_INT, local_buffer, *messageSize, MPI_INT, root, MPI_COMM_WORLD);


    if (*remainder > 0 && *my_rank == root) {
        std::copy(readBuffer + ((*size / *unit) - *remainder), readBuffer + (*size / *unit),
                  local_buffer + *messageSize);
        *messageSize += *remainder;
    }
}


void build_mpi_data_type(int* data_1, int* data_2, int root){

    MPI_Datatype custom_type = NULL;
    // Variables
    MPI_Aint data_1_addr, data_2_addr;
    MPI_Get_address(data_1, &data_1_addr);
    MPI_Get_address(data_2, &data_2_addr);


    int array_of_blocklengths[2] = {1, 1};
    MPI_Datatype array_of_types [2] = {MPI_INT, MPI_INT};
    MPI_Aint array_of_displacements [2] = {0, data_2_addr - data_1_addr};
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements, array_of_types, &custom_type);
    MPI_Type_commit(&custom_type);


    MPI_Bcast(data_1,1, custom_type, root, MPI_COMM_WORLD);
    MPI_Type_free(&custom_type);

}


int main(int argc, char *argv[]) {
    //////// Start Clock //////////////
    // Use Chrono for high grade clock
    int time_samples = 5;
    high_resolution_clock::time_point clock_start[time_samples];
    high_resolution_clock::time_point clock_end[time_samples];
    clock(clock_start, &time_samples);

    //////// USER INPUT//////////////
    if (argc != 3) {
        cout << "Error Error" << endl;
        cout << "Please provide: binary data file and interval size" << endl;
        exit(1);
    }
    assert(argc == 3);
    string filePath = argv[1];
    int intervalSize = check_user_number(argv[2]);
    assert(intervalSize > 0);


    //////// MPI  INIT //////////////
    int comm_sz;  // Number of process
    int my_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    //////// MPI  Variables //////////////
    const int root = 0;
    int maxMessage = 1000;
    int* messageSize = &maxMessage;
    int minMessage = 0;
    int local_buffer[maxMessage + comm_sz];
    init_array(local_buffer, maxMessage + comm_sz);
    int local_data[intervalSize];
    init_array(local_data, intervalSize);
    int local_min = numeric_limits<int>::max();
    int local_max = numeric_limits<int>::min();
    int num_iterations;




    //////// Root Variables //////////////
    int unit = sizeof(int);
    int bufferSize = unit * maxMessage * comm_sz;  // read in chuncks of file
    int size = bufferSize;
    int bucketSize = 0;
    int fileLength = 0;
    int min = local_min;
    int max = local_max;
    int readBuffer[bufferSize];
    init_array(readBuffer, bufferSize);
    int data[intervalSize];
    init_array(data, intervalSize);
    int seek = 0;
    int remainder = 0;





    //////// GET File Size //////////////
    if (my_rank == root) {
        cout << "Starting Program: Assignment 2 MPI" << endl;
        //////// OPEN FILE //////////////
        ifstream fileInput;
        fileInput.open(filePath, ios::binary);
        if (fileInput.is_open()) {
            fileInput.seekg(0, ios::end);
            fileLength = fileInput.tellg();
            num_iterations = (int) ceil((double) fileLength / bufferSize);

            if (fileLength < bufferSize) {
                minMessage = (fileLength / unit) / comm_sz;
            } else{
                minMessage = (( fileLength -  (num_iterations - 1) * bufferSize)) / (unit * comm_sz);
                if (minMessage <= 0 )
                    minMessage = maxMessage;
            }


            fileInput.close();
        } else {
            cout << "Can Not open file..." << endl;
            exit(1);
        }
    }

    // Let Process know how many rounds to expect
//    MPI_Bcast(&num_iterations, 1, MPI_INT, root, MPI_COMM_WORLD);

//    build_mpi_data_type(&num_iterations, messageSize, root);
    build_mpi_data_type(&num_iterations, &minMessage, root);


    int last = num_iterations -1;
    //////// FIND MIN MAX //////////////
    for (int i = 0; i < num_iterations; i++) {
        if(i == last) {
            messageSize = &minMessage;
//            cout << " Message size Min Max: " << *messageSize << " :from processor: " << my_rank << endl;
        }

        if (my_rank == 0) {
            readFile(&filePath, readBuffer, messageSize, &fileLength, &comm_sz, &seek, &size, &unit, &remainder);
        }

        send_data(readBuffer, local_buffer, messageSize, root, &my_rank, &size, &unit, &remainder);

        //Calculate
        for (int i = 0; i < *messageSize; i++) {
            if (local_min > local_buffer[i])
                local_min = local_buffer[i];
            if (local_max < local_buffer[i])
                local_max = local_buffer[i];
        }
    }

    //////// Rest Values and Send Min and Max //////////////
    messageSize = &maxMessage;
    seek = 0;
    remainder = 0;
    size = bufferSize;


//    MPI_Barrier(MPI_COMM_WORLD);
//    int time_samples = 5;
//    high_resolution_clock::time_point clock_start[time_samples];
//    high_resolution_clock::time_point clock_end[time_samples];
//    clock(clock_start, &time_samples);

    MPI_Allreduce(&local_min, &min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

//    MPI_Reduce(&local_min, &min, 1, MPI_INT, MPI_MIN, root, MPI_COMM_WORLD);
//    MPI_Reduce(&local_max, &max, 1, MPI_INT, MPI_MAX, root, MPI_COMM_WORLD);
//    build_mpi_data_type(&min, &max, root);


//    MPI_Barrier(MPI_COMM_WORLD);
//    ////////  END CLOCK //////////////
//    //////// GET TIME //////////////
//    if (my_rank == root) {
//        print(fileLength / unit, min, max, bucketSize, data, my_rank, intervalSize);;
//        clock(clock_end, &time_samples);
//        double total_time = calculate_time(clock_start, clock_end, &time_samples);
//        cout << "AVG Time: " << total_time << " Milli Seconds" << endl;
//    }

    //////// FIND Buckets //////////////
    int range = abs(max - min);
    bucketSize = range / intervalSize;
    bucketSize++;


    for (int i = 0; i < num_iterations; i++) {
        if(i == last) {
            messageSize = &minMessage;
//            cout << " Message size Bucket: " << *messageSize << " :from processor: " << my_rank << endl;
        }

        if (my_rank == 0) {
            readFile(&filePath, readBuffer, messageSize, &fileLength, &comm_sz, &seek, &size, &unit, &remainder);
        }

        send_data(readBuffer, local_buffer, messageSize, root, &my_rank, &size, &unit, &remainder);

        //Calculate
        for (int i = 0; i < *messageSize; i++) {
            local_data[(local_buffer[i] - min) / bucketSize]++;
        }
    }


//    for (int i = 0; i < intervalSize; i++) {
//        MPI_Reduce(&local_data[i], &data[i], intervalSize, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
//    }

    MPI_Reduce(&local_data, &data, intervalSize, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

    MPI_Finalize();

//    ////////  END CLOCK //////////////
//    //////// GET TIME //////////////
    if (my_rank == root) {
        print(fileLength / unit, min, max, bucketSize, data, my_rank, intervalSize);;
        clock(clock_end, &time_samples);
        double total_time = calculate_time(clock_start, clock_end, &time_samples);
        cout << "AVG Time: " << total_time << " Milli Seconds" << endl;
    }
    return 0;
}