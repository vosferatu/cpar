#include <stdio.h>

#define TRUE 1
#define FALSE 0

struct struct_matrix {
    int nrow;
    int ncol;
    int **data;
};

typedef struct struct_matrix Matrix; 

void create_matrix(Matrix *m, int nrow, int ncol);
void populate_matrix(Matrix *m);
//void print_matrix(Matrix *m, char iden);
void shift_matrix_left(Matrix *m, int block_sz, int initial);
void shift_matrix_up(Matrix *m, int block_sz, int initial);
void matrix_product(Matrix *c, Matrix *a, Matrix *b);
int* create_array_as_matrix(int r, int c);
void populate_array_as_matrix(int *arr, int r, int c);
int array_as_matrix_equals(int *a, int *b, int r, int c);
void print_matrix(int *matrix, int n);
void copy_block(int *src, int *dest, int n);
void copy_blocks(int **src, int **dest, int m, int n);
int** create_block(int blocks_pp, int block_sz);
void shift_left(int **buf, int bpp, int bz);
void shift_up(int **buf, int bpp, int bz);
void print_blocks(int **buf, int bpp, int bz);