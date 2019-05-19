#include <iostream>
#include <cstdlib>
#include <chrono>
#include <omp.h>
#include <cmath>

using namespace std;
using namespace std::chrono;

void block_matrix_mul_parallel(float **A, float **B, float **C, int size, int block_size, int num_threads);

int main(int argc, char **argv)
{
	int size = atoi(argv[1]), block_size = 16, num_threads = 0;
	float **A = new float*[size];
	float **B = new float*[size];
	float **C = new float*[size];
	float **CC = new float*[size];
	cout << "Init has begun" << endl;
	for (int i = 0; i < size; i++)
	{
		A[i] = new float[size];
		B[i] = new float[size];
		C[i] = new float[size];
		CC[i] = new float[size];
		for (int j = 0; j < size; j++)
		{
			A[i][j] = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);
			B[i][j] = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);
			C[i][j] = 0.0f;
			CC[i][j] = 0.0f;
		}
	}


	for(int thread_count = 2; thread_count <= 64; thread_count += num_threads,num_threads=thread_count){
		//parallel block
		cout << "----------------------------------\nComputation has begun(parallel block)" << endl;
		cout << "Size: " << size << "     Num_threads: " << num_threads << endl;
    	high_resolution_clock::time_point t1 = high_resolution_clock::now();

		block_matrix_mul_parallel(A, B, CC, size, block_size, num_threads);

    	high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		cout << "Time elapsed: " << static_cast <float> (duration) / 1000000.0f << endl << endl;
	}

	return 0;
}

void block_matrix_mul_parallel(float **A, float **B, float **C, int size, int block_size, int num_threads)
{
	int i = 0, j = 0, k = 0, jj = 0, kk = 0;
	float tmp;
	int chunk = 1;
	int tid;
	omp_set_dynamic(0);
	omp_set_num_threads(num_threads);

#pragma omp parallel shared(A, B, C, size, chunk) private(i, j, k, jj, kk, tid, tmp)
	{
		tid = omp_get_thread_num();
		if (tid == 0)
		{
			cout << "Number of threads: " << omp_get_num_threads() << endl;
		}
		#pragma omp for schedule (static, chunk)
		for (jj = 0; jj < size; jj += block_size)
		{
			for (kk = 0; kk < size; kk += block_size)
			{
				for (i = 0; i < size; i++)
				{
					for (j = jj; j < ((jj + block_size) > size ? size : (jj + block_size)); j++)
					{
						tmp = 0.0f;
						for (k = kk; k < ((kk + block_size) > size ? size : (kk + block_size)); k++)
						{
							tmp += A[i][k] * B[k][j];
						}
						C[i][j] += tmp;
					}
				}
			}
		}
	}
}
