#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "symnmf.h"

/*region VARS*/
int dimensionOfVector = 0;
int numberOfVectors = 0;
int K = 0;
/*endregion VARS*/

/*region PUT_INPUT_IN_ARRAY*/
short didWeGetFirstLine;
Node *createNode(char data)
{
    Node *newNode = (Node *)malloc(sizeof(Node));
    if (newNode == NULL)
    {
        printf("An Error Has Occurred 1");
        exit(1);
    }
    newNode->data = data;
    newNode->prev = NULL;
    newNode->next = NULL;
    return newNode;
}

/* Function to initialize the linked list */
void initializeList(LinkedList *list)
{
    list->head = NULL;
    list->tail = NULL;
    list->length = 0;
}

/* Function to insert a node at the end of the list */
void insert(LinkedList *list, char data)
{
    Node *newNode = createNode(data);
    if (list->head == NULL)
    {
        list->head = newNode;
        list->tail = newNode;
    }
    else
    {
        list->tail->next = newNode;
        newNode->prev = list->tail;
        list->tail = newNode;
    }
    list->length++;
}

/* Function to free the memory allocated for the linked list */
void freeList(LinkedList *list)
{
    Node *current = list->head;
    while (current != NULL)
    {
        Node *temp = current;
        current = current->next;
        free(temp);
    }
    list->head = NULL;
    list->tail = NULL;
    list->length = 0;
}

/* Function to print the linked list */
void printList(LinkedList *list)
{
    Node *current = list->head;
    printf("Normal List: ");
    while (current != NULL)
    {
        printf("%c ", current->data);
        current = current->next;
    }
    printf("\n");

    current = list->tail;
}

/*function to read the file and put the input in an array*/
ArrayInfo read_file_to_array(char *filename)
{
    int wordlength;
    int i, j, k;
    double **array;
    char *word;
    char ch;
    Node *current;
    Node *wordStart;
    FILE *file = fopen(filename, "r");
    ArrayInfo result;
    LinkedList allChars;
    initializeList(&allChars);
    dimensionOfVector = 1;
    numberOfVectors = 0;
    didWeGetFirstLine = 0;

    if (file == NULL)
    {
        printf("An Error Has Occurred 2");
        exit(1);
    }

    /*read the file and put the input in an array*/
    while ((ch = fgetc(file)) != EOF)
    {
        if (ch == ',' && didWeGetFirstLine == 0)
        {
            dimensionOfVector++;
        }
        if (ch == '\n')
        {
            insert(&allChars, ',');
            numberOfVectors++;
            didWeGetFirstLine = 1;
        }
        else
        {
            insert(&allChars, ch);
        }
    }
    fclose(file); /*Close the file*/
    array = malloc(numberOfVectors * sizeof(double *));
    if (array == NULL)
    {
        printf("An Error Has Occurred 3");
        exit(1);
    }
    current = allChars.head;
    wordStart = current;
    wordlength = 0;

    for (i = 0; i < numberOfVectors; i++)
    {
        array[i] = malloc(dimensionOfVector * sizeof(double));
        if (array[i] == NULL)
        {
            printf("An Error Has Occurred 4");
            exit(1);
        }

        for (j = 0; j < dimensionOfVector; j++)
        {
            wordlength = 0;
            wordStart = current;
            while (current->data != ',')
            {
                current = current->next;
                wordlength++;
            }
            current = current->next;

            word = malloc((wordlength + 1) * sizeof(char));
            if (word == NULL)
            {
                for (k = 0; k < i; k++)
                {
                    free(array[k]);
                }
                printf("An Error Has Occurred 5");
                exit(1);
            }
            for (k = 0; k < wordlength; k++)
            {
                word[k] = wordStart->data;
                wordStart = wordStart->next;
            }
            word[wordlength] = '\0'; /*Set the null terminator*/
            array[i][j] = atof(word);
            free(word);
        }
    }

    freeList(&allChars);
    result.array = array;
    result.rows = numberOfVectors;
    result.cols = dimensionOfVector;
    return result;
}
/*endregion PUT_INPUT_IN_ARRAY*/

/*region goals functions*/
/*Function to calculate and output the similarity matrix*/
double **sym(double **X, int rows, int cols)
{
<<<<<<< HEAD
    double **A = allocateMatrix(rows, rows);
    for (int i = 0; i < rows; i++)
=======

    double **A = malloc(rows * sizeof(double *));
    int i;
    int j;
    for (i = 0; i < rows; i++)
    {
        A[i] = malloc(rows * sizeof(double));
    }

    for (i = 0; i < rows; i++)
>>>>>>> ea54a255ef7fb9271ea136fc4fbe7518e3c5d1b0
    {

        for (j = 0; j < rows; j++)
        {
            if (i != j)
            {
                double squaredDistance = 0.0;
                int d;
                for (d = 0; d < cols; d++)
                {
                    double diff = X[i][d] - X[j][d];
                    squaredDistance += diff * diff;
                }
                A[i][j] = exp(-(squaredDistance / 2));
            }
            else
            {
                A[i][j] = 0.0;
            }
        }
    }

    return A;
}

/*Function to calculate and output the diagonal degree matrix*/

double **ddg(double **A, int n)
{
    double **D = malloc(n * sizeof(double *));
    int i, j;
    for (i = 0; i < n; i++)
    {
        D[i] = malloc(n * sizeof(double));
        /*Initialize the row with zeros*/

        for (j = 0; j < n; j++)
        {
            D[i][j] = 0.0;
        }
    }

    /*Calculate the diagonal degree values*/

    for (i = 0; i < n; i++)
    {
        double degree = 0.0;
        for (j = 0; j < n; j++)
        {
            degree += A[i][j];
        }
        D[i][i] = degree;
    }

    return D;
}

/*Allocate memory for a matrix*/
double **allocateMatrix(int rows, int cols)
{
    double **matrix = (double **)malloc(rows * sizeof(double *));
    int i;
    for (i = 0; i < rows; ++i)
    {
        matrix[i] = (double *)malloc(cols * sizeof(double));
    }
    return matrix;
}

/*Calculate the D^-1/2 matrix*/
double **computeDHalfInverse(double **D, int n)
{
    double **D_half_inv = allocateMatrix(n, n);
    int i, j;
    for (i = 0; i < n; ++i)
    {

        for (j = 0; j < n; ++j)
        {
            D_half_inv[i][j] = (i == j) ? 1 / sqrt(D[i][j]) : 0;
        }
    }
    return D_half_inv;
}
/*Multiply two matrices*/
double **matrixMultiply(double **A, double **B, int n)
{
    double **result = allocateMatrix(n, n);
    int i, j, k;
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            result[i][j] = 0;
            for (k = 0; k < n; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}
/*Function to calculate W = D^-1/2 * A * D^-1/2*/
double **norm(double **A, int n)
{
    double **D = ddg(A, n);
    double **D_half_inv = computeDHalfInverse(D, n);
    double **temp = matrixMultiply(D_half_inv, A, n);
    double **W = matrixMultiply(temp, D_half_inv, n);

    /*Free allocated memory for intermediate matrices*/
    int i;
    for (i = 0; i < n; ++i)
    {
        free(D_half_inv[i]);
        free(temp[i]);
    }
    free(D_half_inv);
    free(temp);

    return W;
}

double **symnmf(double **H, double **W, int n, int k)
{
    /*Allocate space for the new H matrix*/
    double **new_H = allocateMatrix(n, k);

    /*Step 1.4.2: Update H*/
    int i, iter, j, l, m;
    double sum_inner;
    for (iter = 0; iter < MAX_ITER; iter++)
    {
        double diffNorm = 0.0;

        /*Update new_H based on current H*/

        for (i = 0; i < n; i++)
        {

            for (j = 0; j < k; j++)
            {
                double numerator = 0.0;
                double denominator = 0.0;

                /*Calculate (WH)ij*/

                for (l = 0; l < n; l++)
                {
                    numerator += W[i][l] * H[l][j];
                }

                /*Calculate (H(H^T)H)ij*/

                for (l = 0; l < n; l++)
                {
                    sum_inner = 0.0;
                    for (m = 0; m < k; m++)
                    {
                        sum_inner += H[i][m] * H[l][m];
                    }
                    denominator += sum_inner * H[l][j];
                }

                new_H[i][j] = H[i][j] * (1 - BETA + BETA * (numerator / (denominator)));
            }
        }

        /*Calculate the norm of the difference between H and new_H*/
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < k; j++)
            {
                diffNorm += pow(new_H[i][j] - H[i][j], 2);
                H[i][j] = new_H[i][j]; /*Update H matrix*/
            }
        }

        /*Check for convergence*/
        if (diffNorm < EPSILON)
        {
            break;
        }
    }

    for (i = 0; i < n; i++)
    {
        free(new_H[i]);
    }
    free(new_H);

    return H;
}

/*endregion goals functions*/

/*region MAIN*/
int main(int argc, char *argv[])
{

    /* DECLARATION ON VARIABLES */
    double **datapoints;
    char *mode = NULL, *input_file_name = NULL;
    int number_datapoints;
    double **outputmatrix;
    double **tempmatrix = NULL;
    int i, j;
    ArrayInfo datapointsstruct;

    if (argc >= 3)
    {
        mode = argv[1];
        input_file_name = argv[2];
    }
    datapointsstruct = read_file_to_array(input_file_name);
    datapoints = datapointsstruct.array;

    number_datapoints = datapointsstruct.rows;

    /*if mode is equal to sym*/
    if (strcmp(mode, "sym") == 0)
    {
        outputmatrix = sym(datapoints, number_datapoints, dimensionOfVector);
    }
    else if (strcmp(mode, "ddg") == 0)
    {
        tempmatrix = sym(datapoints, number_datapoints, dimensionOfVector);
        outputmatrix = ddg(tempmatrix, number_datapoints);
    }
    else if (strcmp(mode, "norm") == 0)
    {
        tempmatrix = sym(datapoints, number_datapoints, dimensionOfVector);
        outputmatrix = norm(tempmatrix, number_datapoints);
    }
    else
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    /*print outputmatrix*/

    for (i = 0; i < number_datapoints; i++)
    {
        for (j = 0; j < number_datapoints; j++)
        {
            if (j != 0)
                printf(",");
            printf("%.4f", outputmatrix[i][j]);
        }
        printf("\n");
    }
    return 0;
}

/*endregion MAIN*/