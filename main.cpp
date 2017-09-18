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

double calculate_time(high_resolution_clock::time_point *start, high_resolution_clock::time_point *end, int *time_samples) {
    // Average time and convert to Micro Sec; 1 sec = 1,000,000 micro sec
    double total = 0;
    for (int i = 0; i < *time_samples; i++) {
        chrono::duration<double, std::milli> diff = end[i] - start[i];  // Time in Micro Sec
        total += diff.count();
    }
    return total / *time_samples;
}

void print(int sample_size, int min, int max, int bucketSize, int* data, int myrank, int intervalSize) {

    cout << "Report from processor:  " << myrank << endl;
    cout << "Sample Size: " << sample_size << endl;
    cout << "Min Value: " <<  min << endl;
    cout << "Max Value: " << max << endl;
    cout << "bucket Range: " << bucketSize - 1 << endl;

    for(int i = 0; i< intervalSize; i++){
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


void min_max_reduce(int *min, int *max, int *local_buffer, int &maxMessage, int root) {
        // todo make a function
}

void init_array(int* a, int arr_size){
    new int[arr_size];
    for(int i = 0; i< arr_size; i++)
        a[i] = 0;
}



void readFile(string* filePath, int* readBuffer, int* messageSize, int* fileLength, int* comm_sz, int* seek, int* size, int* unit){
    //////// OPEN FILE //////////////
    ifstream fileInput;
    fileInput.open(*filePath, ios::binary);
    if (fileInput.is_open()) {
        fileInput.seekg(*seek);

        //  Check if buffer is less then remainder of file
        if (*fileLength - *seek < *size) {
            *size = *fileLength - *seek;
            *messageSize = (*size/ *unit ) / *comm_sz;


            cout << "Changing Message Size"  << endl;
            cout << "File Length: " << *fileLength << endl;
            cout << "Size: " << *size << endl;
            cout << "Seek: " << *seek << endl;
            cout << "Message Size: " << *messageSize << endl;
        }

        fileInput.read((char *) &readBuffer[0], *size);
        *seek = *size;
        for(int i = 0; i < *size/ *unit; i++){
            cout << i <<  " Value: " << readBuffer[i] << endl;

        }

//        *seek = fileInput.tellg();
        fileInput.close();
    } else {
        cout << "Can Not open file..." << endl;
        exit(1);
    }
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



    //////// MPI  Variables //////////////
    const int root = 0;
    int maxMessage = 1000;
    int messageSize = maxMessage;
    int local_buffer[maxMessage];
    init_array(local_buffer, maxMessage);
    int local_data[intervalSize];
    init_array(local_data, intervalSize);
    int local_min = numeric_limits<int>::max();
    int local_max = numeric_limits<int>::min();
    int num_iterations;
    int remainder = 0;


    //////// MPI  INIT //////////////
    int comm_sz;  // Number of process
    int my_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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
            cout << "Num of Iterations: " << num_iterations << endl;
            fileInput.close();
        } else {
            cout << "Can Not open file..." << endl;
            exit(1);
        }
    }

    // Let Process know how many rounds
    MPI_Bcast(&num_iterations, 1, MPI_INT, root, MPI_COMM_WORLD );


    //////// FIND MIN MAX //////////////
    for (int i = 0; i < num_iterations; i++) {
        if (my_rank == 0) {
            readFile(&filePath, readBuffer, &messageSize, &fileLength, &comm_sz, &seek, &size, &unit);
            remainder = (size/ unit) - (messageSize * comm_sz);
            cout << "Remainder Size: " << remainder << " From Processor: " << my_rank << endl;
        }


        MPI_Bcast(&messageSize, 1, MPI_INT, root, MPI_COMM_WORLD );
        MPI_Scatter(&readBuffer, messageSize, MPI_INT, &local_buffer, messageSize, MPI_INT, root, MPI_COMM_WORLD);


        if(remainder > 0 && my_rank == root){
            std::copy(readBuffer + ((size/ unit )- remainder), readBuffer + (size/unit), local_buffer + messageSize);
            messageSize += remainder;
        }






//        // Message length can vary: Update then send array
//        MPI_Bcast(&messageSize, 1, MPI_INT, root, MPI_COMM_WORLD );
//        MPI_Bcast(&remainder, 1, MPI_INT, root, MPI_COMM_WORLD );
//
//
//        if(remainder <= 0){
//            MPI_Scatter(&readBuffer, messageSize, MPI_INT, &local_buffer, messageSize, MPI_INT, root, MPI_COMM_WORLD);
//        } else{
////            int remainder = (messageSize * comm_sz) - (size/ unit);
//            if(my_rank == root){
//                int offset = 0;
//                for(int i = 1; i < comm_sz; i++) {
//                    if (remainder > 0) {
//                        offset = 1;
//                        remainder--;
//                    } else {
//                        offset = 0;
//                    }
//                    int new_message_size = messageSize + offset;
//                    MPI_Send(&readBuffer, new_message_size, MPI_INT, i, 0, MPI_COMM_WORLD);
//
//                }
//                std::copy(readBuffer + size - messageSize, readBuffer + size, local_buffer);
//
//
//            } else{
//                if(my_rank <= remainder)
//                    messageSize++;
//                MPI_Recv(&local_buffer, messageSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//            }
//
//
//        }



        cout << "Message Size: " << messageSize << " From Processor: " << my_rank << endl;


        //Calculate
        for (int i = 0; i < messageSize; i++) {
            if (my_rank == root)
                 cout << "Value: " << local_buffer[i] << " From Processor: " << my_rank << endl;
            if (local_min > local_buffer[i])
                local_min = local_buffer[i];
            if (local_max < local_buffer[i])
                local_max = local_buffer[i];
        }
    }

    //////// Rest Values and Send Min and Max //////////////
    messageSize = maxMessage;
    seek = 0;
    remainder =0;
    size = bufferSize;
    MPI_Allreduce(&local_min, &min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);



    //////// FIND Buckets //////////////
    int range = abs(max - min);
    bucketSize =  range / intervalSize;
    bucketSize++;

//    for(int i = 0; i< intervalSize; i++){
//        cout << "Value: " << local_data[i] << " From Process: " << my_rank <<endl;
//    }



    for (int i = 0; i < 1; i++) {
        if (my_rank == 0) {
            readFile(&filePath, readBuffer, &messageSize, &fileLength, &comm_sz, &seek, &size, &unit);
            remainder = (size/ unit) - (messageSize * comm_sz);
        }



//        //////////  Temp
//
//        // Message length can vary: Update then send array
//        MPI_Bcast(&messageSize, 1, MPI_INT, root, MPI_COMM_WORLD );
//        MPI_Bcast(&remainder, 1, MPI_INT, root, MPI_COMM_WORLD );
//
//
//        if(remainder <= 0){
//            MPI_Scatter(&readBuffer, messageSize, MPI_INT, &local_buffer, messageSize, MPI_INT, root, MPI_COMM_WORLD);
//        } else{
////            int remainder = (messageSize * comm_sz) - (size/ unit);
//            if(my_rank == root){
//                int offset = 0;
//                for(int i = 1; i < comm_sz; i++) {
//                    if (remainder > 0) {
//                        offset = 1;
//                        remainder--;
//                    } else {
//                        offset = 0;
//                    }
//                    int new_message_size = messageSize + offset;
//                    MPI_Send(&readBuffer, new_message_size, MPI_INT, i, 0, MPI_COMM_WORLD);
//
//                }
//                std::copy(readBuffer + size - messageSize, readBuffer + size, local_buffer);
//
//
//            } else{
//                if(my_rank <= remainder)
//                    messageSize++;
//                MPI_Recv(&local_buffer, messageSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//            }
//
//
//        }
//
//
//
//
//        /////////  End Temp














        // Message length can vary: Update then send array
        MPI_Bcast(&messageSize, 1, MPI_INT, root, MPI_COMM_WORLD );
//        cout << "Message Size: " << messageSize << " From Process: " << my_rank <<endl;
        MPI_Scatter(&readBuffer, messageSize, MPI_INT, &local_buffer, messageSize, MPI_INT, root, MPI_COMM_WORLD);

        if(remainder > 0 && my_rank == root){
            std::copy(readBuffer + ((size/ unit )- remainder), readBuffer + (size/unit), local_buffer + messageSize);
            messageSize += remainder;
        }

        //Calculate





        for (int i = 0; i < messageSize; i++) {
            local_data[(local_buffer[i] - min) / bucketSize]++;
        }
    }
//    for(int i = 0; i< intervalSize; i++){
//        cout << "Value: " << local_data[i] << " From Process: " << my_rank <<endl;
//    }

    for(int i = 0; i< intervalSize; i++){
        MPI_Reduce(&local_data[i], &data[i], 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    }


    MPI_Finalize();



    ////////  END CLOCK //////////////
    //////// GET TIME //////////////
    if(my_rank == root){
        print(fileLength / unit, min, max, bucketSize, data, my_rank, intervalSize);;
        clock(clock_end, &time_samples);
        double total_time = calculate_time(clock_start, clock_end, &time_samples);
        cout << "AVG Time: " << total_time << " Milli Seconds" << endl;
    }


    return 0;
}