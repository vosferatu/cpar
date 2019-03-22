#include <stdio.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <chrono>
#include <algorithm>

using namespace std;

#define SYSTEMTIME clock_t

double * gen_matrix(int size, bool grad){

	double * res = (double *)malloc((size * size) * sizeof(double));

	for(int i = 0; i < size; i++)
		for(int j = 0; j  < size; j++)
			if(grad)
				res[i*size + j] = j+i*size;
			else
				res[i*size + j] = 1;

	return res;

}


double *  ijk_implementation(double * A, double * B, double * C, int size){

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                C[i*size + j] += A[i*size + k] * B[k*size + j];
            }
        }
    }

    return C;
} 

double * ikj_implementation(double * A, double * B, double * C, int size){

    for (int i = 0; i < size; i++) {
        for (int k = 0; k < size; k++) {
            for (int j = 0; j < size; j++) {
                C[i*size + j] += A[i*size + k] * B[k*size + j];
            }
        }
    }
    return C;
} 


double * matrix_block(double * A, double * B, double * C, int size){

    int i,j,k,x,y,z;

// tiled (parallel) matrix multiplication
// from /sys/devices/system/cpu/cpu0/cache/index0
// cat c
    
    const char* env_p = std::getenv("L1DSIZE");

    int incr = atoi(env_p);

    for (i = 0; i < size; i += incr) {
         for (j = 0; j < size; j += incr) {
             C[i*size+j] = 0.0;
             for (k = 0; k < size; k += incr){
             	for (z = k; z < std::min( k + incr, size ); z++) {
                	 for (x = i; x < std::min( i + incr, size ); x++) {
                    	 for (y = j; y < std::min( j + incr, size ); y++) {
                 
                             C[ x * size + y ] +=  A[ x * size + z ] * B[ z * size  + y  ];

                         }
                     } 
                 }
             }
         }
     }


     return C;


}




bool equal_matrix(double * A, double* B, int size){
	for(int i=0; i<size ; i++)
		for(int j=0; j<size; j++)
			if(A[i*size + j] != B[i*size + j])
				return false;

	return true;
}

void print_matrix(double * A, int size){
	cout << "\n Printing Matrix \n";

	for(int i=0; i<size ; i++){
		for(int j=0; j<size; j++){
			cout << ", " << A[i*size + j] << " ";
		}
		cout << " \n ";
	}
	cout << " \n";
}


bool test_implementation(double * (*implementation)(double * A, double * B, double * C, int size), int size){

	double * (*correct_implementation)(double * A1, double * B1, double * C1, int size1) = &ikj_implementation;
	
	cout << "Runnning tests for implementation. Seeing if correctly implemented.  \n.";

	double * test_matrix = gen_matrix(size, true);

	double * correct_results = (double *) malloc((size*size)* sizeof(double));
	
	double * test_results = (double *) malloc((size*size)* sizeof(double));

	correct_implementation(test_matrix, test_matrix, correct_results, size);	


	implementation(test_matrix, test_matrix, test_results, size);


	if(memcmp(correct_results, test_results, size*size*sizeof(double)) == 0){
		cout << "OK, matrix equal \n";
		return 0;
	}
	else{
		cout << "SOMETHING Wrong, matrix not equal \n";		
		return -1;
	}

}

void time_implementation(double * (*implementation)(double * A, double * B, double * C, int size), int size, int run_times){
	
	cout << "Timing for implementation.  \n";

	double * test_matrix = gen_matrix(size, true);

	double * test_results = (double *) malloc((size*size)* sizeof(double));

	auto start = std::chrono::steady_clock::now();

	for(int i = 0; i < run_times; i++ ){
		memset(test_results, 0, (size*size)* sizeof(double));
		implementation(test_matrix, test_matrix, test_results, size);

	}

	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now() - start).count();

	elapsed = elapsed*1.0/run_times;

	cout << "Implementation took: " << elapsed << "ms to run on matrix of size: " << size << "x" << size << "\n";		
}	

 bool ispowerof2(int x) {
   return x && !(x & (x - 1));
 }


int main(int argc, char const *argv[])
{

	int size;
	bool testing = false;
	double * (*imp)(double * A1, double * B1, double * C1, int size1) = &matrix_block;

	if(argc < 3 || argc > 4){
		cout << "Incorrect number of arguments, specify size and if test or not. \n";
		return -1;
	} else {


		if(strcmp("block", argv[1]) == 0)
			imp = &matrix_block;

		if(strcmp("ijk", argv[1]) == 0)
			imp = &ijk_implementation;

		if(strcmp("ikj", argv[1]) == 0)
			imp = &ikj_implementation;

		size = atoi(argv[2]);

		if(!(size > 0 && size < 16384 )){
			cout << "Size needs to be between 0 and 16384 and power of 2. \n";
			return -1;
		}

		if(argc == 4){
			testing = true;
		}
	}


	if(testing)
		time_implementation(imp, size, 3);
	else
		time_implementation(imp, size, 3);
	return 0;
}
