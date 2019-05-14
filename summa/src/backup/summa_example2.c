#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "c_timer.h"
#include "mpi.h"

int main(int argc, char* argv[]) {
	srand((unsigned int)time(NULL));
	int LDA=0,
		LDB=0,
		LDC=0,
		N=0,
		M=0,
		P=0,
		NTROW=0,
		NTCOL=0,
		dimension=0,
		SR=0,
		SC=0,
		rank=0,
		nproc=0,
		bl_sz_n=0,
		bl_sz_m=0,
		bl_sz_p=0;
	double inizio, fine, tempo, tempo_totale;

	void nullifyMatrix(int LD, float M[][LD], int N),
		printMatrix(int LD, float M[][LD], int N),
		randomMatrix(int LD, float M[][LD], int N),
		summa(int LDA, int LDB, int LDC, float A[][LDA], float B[][LDB], float C[][LDC], int Nb, int Mb, int Pb, int SR, int SC, int NTROW, int NTCOL, MPI_Comm comm);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	if(argc == 4) {
		dimension = atoi(argv[1]);
		NTROW = atoi(argv[2]);
		NTCOL = atoi(argv[3]);
	} else {
		if(rank == 0)
			perror("./summa [LD] [NTROW] [NTCOL]");
		MPI_Finalize();
		return 0;
	}

	LDA = LDB = LDC = N = M = P = dimension;
	SR = SC = sqrt(nproc);
	bl_sz_n = N / SR;
	bl_sz_m = M / SC;
	bl_sz_p = P / SC;

	LDA = bl_sz_n;
	LDB = bl_sz_m;
	LDC = bl_sz_p;

	/*
	** IPOTIZZO CHE OGNI PROCESSORE SI CREI LA SUA SOTTOMATRICE CON VALORI RANDOM.
	** NELLA REALTA' OGNUNO RICEVE LA SUA SOTTOMATRICE IN INPUT
	*/
	float *A_loc = (float *) malloc(bl_sz_n * bl_sz_m * sizeof(float));
	float *B_loc = (float *) malloc(bl_sz_m * bl_sz_p * sizeof(float));
	float *C_loc = (float *) malloc(bl_sz_n * bl_sz_p * sizeof(float));

	randomMatrix(bl_sz_n, (float (*)[LDA]) A_loc, bl_sz_n);
	randomMatrix(bl_sz_m, (float (*)[LDB]) B_loc, bl_sz_m);
	nullifyMatrix(bl_sz_p, (float (*)[LDC]) C_loc, bl_sz_p);

	inizio = get_cur_time();
	summa(LDA, LDB, LDC, (float (*)[LDA]) A_loc, (float (*)[LDB]) B_loc, (float (*)[LDC]) C_loc,
		bl_sz_n, bl_sz_m, bl_sz_p, SR, SC, NTROW, NTCOL, MPI_COMM_WORLD);
	fine = get_cur_time();
	tempo = fine - inizio;
	
	MPI_Reduce(&tempo, &tempo_totale, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);	
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0) printf("Total elapsed time: %lfs\n", tempo_totale);

	MPI_Finalize();
	return 0;
}

void printMatrix(int LD, float M[][LD], int N) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			printf("[%2.1f] ", M[i][j]);
		}
		printf("\n");
	}
}

void randomMatrix(int LD, float M[][LD], int N) {
	int i, j, m=0;
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			m = rand() % 100 + 1;
			M[i][j]=m;
		}
	}
}

void nullifyMatrix(int LD, float M[][LD], int N) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			M[i][j]=0;
		}
	}
}

void matmatikj(int LDA, int LDB, int LDC, float A[][LDA], float B[][LDB], float C[][LDC], int N, int M, int P) {
	int i=0, j=0, k=0;
	for(i=0; i<N; i++)
		for(k=0; k<M; k++)
			for(j=0; j<P; j++)
				C[i][j] += A[i][k] * B[k][j];
}

void matmatthread(int LDA, int LDB, int LDC, float A[][LDA], float B [][LDB], float C[][LDC], int N, int M, int P, int NTROW, int NTCOL) {
	int n_threads = NTROW * NTCOL, i=0, j=0;
	pthread_t **threads = (pthread_t **) malloc(NTROW * sizeof(pthread_t *));
	void *thread_function(void *arg);

	typedef struct {
		int row, column, LDA, LDB, LDC, N, M, P;
		float *A, *B, *C;
	} parameters;

	for(i=0; i < NTROW; i++)
		threads[i] = (pthread_t *) malloc(NTCOL * sizeof(pthread_t));

	N /= NTROW;
	P /= NTCOL;

	parameters *param;
	for(i=0; i < NTROW; i++) {
		for(j=0; j < NTCOL; j++) {
			param = (parameters *) malloc(sizeof(parameters));
			param->row = i;
			param->column = j;
			param->LDA = LDA;
			param->LDB = LDB;
			param->LDC = LDC;
			param->N = N;
			param->M = M;
			param->P = P;
			param->A = *A;
			param->B = *B;
			param->C = *C;
			pthread_create(&threads[i][j], NULL, thread_function, param);
		}
	}
	for(i=0; i < NTROW; i++) {
		for(j=0; j < NTCOL; j++) {
			pthread_join(threads[i][j], NULL);
		}
	}
}

void *thread_function(void *arg) {
	void matmatikj(int LDA, int LDB, int LDC, float A[][LDA], float B[][LDB], float C[][LDC], int N, int M, int P);

	typedef struct {
		int row, column, LDA, LDB, LDC, N, M, P;
		float *A, *B, *C;
	} parameters;

	parameters *param = (parameters *) arg;
	int row = param->row,
		column = param->column,
		LDA = param->LDA,
		LDB = param->LDB,
		LDC = param->LDC,
		N = param->N,
		M = param->M,
		P = param->P;
	float *A = param->A,
		*B = param->B,
		*C = param->C;

	float *myA = (float *) A + row * N * LDA, *myB = (float *) B + column * P, *myC = (float *) C + row * N * LDC + column * P;
	matmatikj(LDA, LDB, LDC, (float (*)[LDA]) myA, (float (*)[LDB]) myB, (float (*)[LDC]) myC, N, M, P);
}

void summa(int LDA, int LDB, int LDC, float A[][LDA], float B[][LDB], float C[][LDC], int Nb, int Mb, int Pb, int SR, int SC, int NTROW, int NTCOL, MPI_Comm comm) {
	int rank=0;
	MPI_Comm_rank(comm, &rank);

	int i=0, k=0, row_color = rank / SR, col_color = rank % SC + SC, bl_sz_A = Nb * Mb, bl_sz_B = Mb * Pb;

	void nullifyMatrix(int LD, float M[][LD], int N),
		matmatthread(int LDA, int LDB, int LDC, float A[][LDA], float B[][LDB], float C[][LDC], int N, int M, int P, int NTROW, int NTCOL);

	float *A_temp = (float *) malloc(bl_sz_A * sizeof(float));
	float *B_temp = (float *) malloc(bl_sz_B * sizeof(float));
	nullifyMatrix(LDA, (float (*)[LDA]) A_temp, Nb);
	nullifyMatrix(LDA, (float (*)[LDA]) B_temp, Mb);

	MPI_Comm row_comm, col_comm;
	MPI_Comm_split(comm, row_color, rank, &row_comm);
	MPI_Comm_split(comm, col_color, rank, &col_comm);

	for (k = 0; k < SR; k++) {
		if (col_color == k + SC) memcpy(A_temp, A, bl_sz_A * sizeof(float));
		if (row_color == k) memcpy(B_temp, B, bl_sz_B * sizeof(float));

		MPI_Bcast(A_temp, bl_sz_A, MPI_FLOAT, k, row_comm);
		MPI_Bcast(B_temp, bl_sz_B, MPI_FLOAT, k, col_comm);
		matmatthread(LDA, LDB, LDC, (float (*)[LDA]) A_temp, (float (*)[LDB]) B_temp, (float (*)[LDC]) C, Nb, Mb, Pb, NTROW, NTCOL);
	}
}
