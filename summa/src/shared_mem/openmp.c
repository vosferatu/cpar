#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <omp.h>
#include <math.h>

typedef double TYPE;
#define MAX_DIM 4096*4096
#define MAX_VAL 10
#define MIN_VAL 1
#define ITERATIONS 1

TYPE** randomSquareMatrix(int dimension);
TYPE** zeroSquareMatrix(int dimension);
void displaySquareMatrix(TYPE** matrix, int dimension);
void convert(TYPE** matrixA, TYPE** matrixB, int dimension);

double parallelMultiply(TYPE** matrixA, TYPE** matrixB, TYPE** matrixC, int dimension, int thread_count);
void parallelMultiplyTest(int dimension, int num_threads);

// 1 Dimensional matrix on stack
TYPE flatA[MAX_DIM];
TYPE flatB[MAX_DIM];

int main(int argc, char* argv[]){
	int num_threads = strtol(argv[1], NULL, 10);

	FILE *fp;
	// Create Parallel Multiply test log
	fp = fopen("ParallelMultiplyTest.txt", "w+");
	fclose(fp);

	int step = 128;
	for(int dimension=128; dimension<=4096; dimension+=step){
		parallelMultiplyTest(dimension, num_threads);
		step = dimension;
	}

	return 0;
}

TYPE** randomSquareMatrix(int dimension){

	TYPE** matrix = malloc(dimension * sizeof(TYPE*));

	for(int i=0; i<dimension; i++){
		matrix[i] = malloc(dimension * sizeof(TYPE));
	}

	//Random seed
	srandom(time(0)+clock()+random());

	#pragma omp parallel for
	for(int i=0; i<dimension; i++){
		for(int j=0; j<dimension; j++){
			matrix[i][j] = rand() % MAX_VAL + MIN_VAL;
		}
	}

	return matrix;
}

TYPE** zeroSquareMatrix(int dimension){

	TYPE** matrix = malloc(dimension * sizeof(TYPE*));

	for(int i=0; i<dimension; i++){
		matrix[i] = malloc(dimension * sizeof(TYPE));
	}

	//Random seed
	srandom(time(0)+clock()+random());
	for(int i=0; i<dimension; i++){
		for(int j=0; j<dimension; j++){
			matrix[i][j] = 0;
		}
	}

	return matrix;
}

void displaySquareMatrix(TYPE** matrix, int dimension){
	for(int i=0; i<dimension; i++){
		for(int j=0; j<dimension; j++){
			printf("%f\t", matrix[i][j]);
		}
		printf("\n");
	}
}

double parallelMultiply(TYPE** matrixA, TYPE** matrixB, TYPE** matrixC, int dimension, int thread_count){

	int i, j, k, iOff, jOff;
	TYPE tot;

	struct timeval t0, t1;
	gettimeofday(&t0, 0);

	convert(matrixA, matrixB, dimension);
	#pragma omp parallel shared(matrixC) private(i, j, k, iOff, jOff, tot) num_threads(thread_count)
	{
		#pragma omp for schedule(static)
		for(i=0; i<dimension; i++){
			iOff = i * dimension;
			for(j=0; j<dimension; j++){
				jOff = j * dimension;
				tot = 0;
				for(k=0; k<dimension; k++){
					tot += flatA[iOff + k] * flatB[jOff + k];
				}
				matrixC[i][j] = tot;
			}
		}
	}

	gettimeofday(&t1, 0);
	double elapsed = (t1.tv_sec-t0.tv_sec) * 1.0f + (t1.tv_usec - t0.tv_usec) / 1000000.0f;

	return elapsed;
}

void convert(TYPE** matrixA, TYPE** matrixB, int dimension){
	#pragma omp parallel for
	for(int i=0; i<dimension; i++){
		for(int j=0; j<dimension; j++){
			flatA[i * dimension + j] = matrixA[i][j];
			flatB[j * dimension + i] = matrixB[i][j];
		}
	}
}

void parallelMultiplyTest(int dimension, int num_threads){
	FILE* fp;
	fp = fopen("ParallelMultiplyTest.txt", "a+");

	// Console write
	printf("----------------------------------\n");
	printf("Test : Parallel Multiply\n");
	printf("----------------------------------\n");
	printf("Dimension : %d\n", dimension);
	printf("..................................\n");

	// File write
	fprintf(fp, "----------------------------------\n");
	fprintf(fp, "Test : Parallel Multiply\n");
	fprintf(fp, "----------------------------------\n");
	fprintf(fp, "Dimension : %d\n", dimension);
	fprintf(fp, "..................................\n");

	double* opmLatency = malloc(ITERATIONS * sizeof(double));
	TYPE** matrixA = randomSquareMatrix(dimension);
	TYPE** matrixB = randomSquareMatrix(dimension);

	// Iterate and measure performance
	for(int i=0; i<ITERATIONS; i++){
		TYPE** matrixResult = zeroSquareMatrix(dimension);
		opmLatency[i] = parallelMultiply(matrixA, matrixB, matrixResult, dimension, num_threads);
		free(matrixResult);

		// Console write
		printf("%d.\t%f\n", i+1, opmLatency[i]);

		// File write
		fprintf(fp, "%d.\t%f\n", i+1, opmLatency[i]);
	}

	// Console write
	printf("\n");
	printf("----------------------------------\n");
	printf("Analyze Measurements              \n");
	printf("----------------------------------\n");

	// File write
	fprintf(fp, "\n");
	fprintf(fp, "----------------------------------\n");
	fprintf(fp, "Analyze Measurements              \n");
	fprintf(fp, "----------------------------------\n");

	double sum = 0.0;
	double sumSquared = 0.0;

	// Statistical analyze
	for(int i=0; i<ITERATIONS; i++){
		sum += opmLatency[i];
		sumSquared += pow(opmLatency[i], 2.0);
	}

	double mean = sum / ITERATIONS;
	double squareMean = sumSquared / ITERATIONS;
	double standardDeviation = sqrt(squareMean - pow(mean, 2.0));

	// Console write
	printf("Mean               : %f\n", mean);
	printf("Standard Deviation : %f\n", standardDeviation);
	printf("----------------------------------\n");

	//File write
	fprintf(fp, "Mean               : %f\n", mean);
	fprintf(fp, "Standard Deviation : %f\n", standardDeviation);
	fprintf(fp, "----------------------------------\n");

	// Releasing memory
	fclose(fp);
	free(opmLatency);
	free(matrixA);
	free(matrixB);
}
