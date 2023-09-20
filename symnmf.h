#ifndef SYM_NMF_H
#define SYM_NMF_H

// Define constants
#define EPSILON 0.0001
#define MAX_ITER 300
#define BETA 0.5

typedef struct Node
{
    char data;
    struct Node *prev;
    struct Node *next;
} Node;

typedef struct
{
    Node *head;
    Node *tail;
    int length;
} LinkedList;

typedef struct
{
    double **array;
    int rows;
    int cols;
} ArrayInfo;

typedef struct
{
    double **array;
    int rows;
    int cols;
} Matrix;

// Function prototypes
ArrayInfo read_file_to_array(char *filename);
double **sym(double **X, int rows, int cols);
double **ddg(double **A, int n);
double **norm(double **A, int n);
double **symnmf(double **H, double **W, int n, int k);
double **multiplyMatrices(double **mat1, int rows1, int cols1, double **mat2, int rows2, int cols2);
double **transposeMatrix(double **mat, int n, int d);
double **calculation(double **H, double **W_H, double **H_HT_H, int H_n, int H_d);
void copyMatrix(double **sourceMatrix, double **destinationMatrix, int rows, int cols);
void freeMatrix(double **matrix, int rows);
double **mallocMatrix(int rows, int cols);
void updateH(double **W, double **old_H, int W_n, int old_H_n, int old_H_d, double **new_H);

#endif /* SYM_NMF_H */
