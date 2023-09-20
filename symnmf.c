#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "symnmf.h"

/*region VARS*/
int dimensionOfVector;
int numberOfVectors;
int K;
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

ArrayInfo read_file_to_array(char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("An Error Has Occurred 2");
        exit(1);
    }

    int wordlength;
    int i, j, k;
    double **array;
    char *word;
    char ch;
    Node *current;
    Node *wordStart;
    LinkedList allChars;
    initializeList(&allChars);
    dimensionOfVector = 1;
    numberOfVectors = 0;
    didWeGetFirstLine = 0;

    while ((ch = fgetc(file)) != EOF) // Changed getchar() to fgetc(file)
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
    fclose(file); // Close the file
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

            /* wordStart; */

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
    ArrayInfo result;
    result.array = array;
    result.rows = numberOfVectors;
    result.cols = dimensionOfVector;
    return result;
}
/*endregion PUT_INPUT_IN_ARRAY*/
// free matrix
void freeMatrix(double **matrix, int rows)
{
    for (int i = 0; i < rows; ++i)
    {
        free(matrix[i]);
    }
    free(matrix);
}

/*region goals functions*/
// Function to calculate and output the similarity matrix
double **sym(double **X, int rows, int cols)
{

    double **A = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        A[i] = malloc(rows * sizeof(double));
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            if (i != j)
            {
                double squaredDistance = 0.0;
                for (int d = 0; d < cols; d++)
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

    return A; // Return the dynamically allocated array
}

// Function to calculate and output the diagonal degree matrix

double **ddg(double **A, int n)
{
    // Allocate memory for the degree matrix D
    double **D = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        D[i] = malloc(n * sizeof(double));
        // Initialize the row with zeros
        for (int j = 0; j < n; j++)
        {
            D[i][j] = 0.0;
        }
    }

    // Calculate the diagonal degree values
    for (int i = 0; i < n; i++)
    {
        double degree = 0.0;
        for (int j = 0; j < n; j++)
        {
            degree += A[i][j];
        }
        D[i][i] = degree;
    }

    return D;
}

// Allocate memory for a matrix
double **allocateMatrix(int rows, int cols)
{
    double **matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; ++i)
    {
        matrix[i] = (double *)malloc(cols * sizeof(double));
    }
    return matrix;
}

// Calculate the D^-1/2 matrix
double **computeDHalfInverse(double **D, int n)
{
    double **D_half_inv = allocateMatrix(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            D_half_inv[i][j] = (i == j) ? 1 / sqrt(D[i][j]) : 0;
        }
    }
    return D_half_inv;
}
// Multiply two matrices
double **matrixMultiply(double **A, double **B, int n)
{
    double **result = allocateMatrix(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i][j] = 0;
            for (int k = 0; k < n; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}
// Function to calculate W = D^-1/2 * A * D^-1/2
double **norm(double **A, int n)
{
    double **D = ddg(A, n);
    double **D_half_inv = computeDHalfInverse(D, n);
    double **temp = matrixMultiply(D_half_inv, A, n);
    double **W = matrixMultiply(temp, D_half_inv, n);

    // Free allocated memory for intermediate matrices
    for (int i = 0; i < n; ++i)
    {
        free(D_half_inv[i]);
        free(temp[i]);
    }
    free(D_half_inv);
    free(temp);

    return W;
}
double euclideanDistance(double **mat1, double **mat2, int rows, int cols)
{
    double sum = 0.0;
    int i;
    int j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            double diff = mat1[i][j] - mat2[i][j];
            sum += diff * diff;
        }
    }
    return sum;
}
double **symnmf(double **H, double **W, int n, int k)
{
    double diff; // arbitrary large number
    int iteration = 0;

    // Allocate space for the new H matrix
    double **new_H = mallocMatrix(n, k);

    // Iterate until convergence or maximum iterations
    while (iteration < MAX_ITER)
    {
        // Update H matrix
        updateH(W, H, n, n, k, new_H);

        // Calculate the difference between new_H and H for convergence check
        diff = euclideanDistance(new_H, H, n, k);
        
        // Otherwise, set H to new_H for next iteration and update previous_diff
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                H[i][j] = new_H[i][j];
            }
        }
        iteration++;

        if (diff < EPSILON)
        {
            printf("iter=%d\n", iteration);
            break;
        }
    }

    freeMatrix(new_H, n); // Free the allocated memory for new_H

    return H; // Return the final H matrix
}

void updateH(double **W, double **old_H, int W_n, int old_H_n, int old_H_d, double **new_H)
{
    double **H_transpose;
    double **W_H;
    double **H_HT;
    double **H_HT_H;
    double **result;

    H_transpose = transposeMatrix(old_H, old_H_n, old_H_d);
    W_H = multiplyMatrices(W, W_n, W_n, old_H, old_H_n, old_H_d);
    H_HT = multiplyMatrices(old_H, old_H_n, old_H_d, H_transpose, old_H_d, old_H_n);
    H_HT_H = multiplyMatrices(H_HT, old_H_n, old_H_n, old_H, old_H_n, old_H_d);
    result = calculation(old_H, W_H, H_HT_H, old_H_n, old_H_d);
    copyMatrix(result, new_H, old_H_n, old_H_d);

    freeMatrix(H_transpose, old_H_d);
    freeMatrix(W_H, W_n);
    freeMatrix(H_HT, old_H_n);
    freeMatrix(H_HT_H, old_H_n);
    freeMatrix(result, old_H_n);
}

// mallocMatrix function
double **mallocMatrix(int rows, int cols)
{
    int i;
    double **matrix;
    matrix = malloc(rows * sizeof(double *));
    if (matrix == NULL)
    {
        printf("An Error Has Occurred 6");
        exit(1);
    }
    for (i = 0; i < rows; i++)
    {
        matrix[i] = malloc(cols * sizeof(double));
        if (matrix[i] == NULL)
        {
            printf("An Error Has Occurred 7");
            exit(1);
        }
    }
    return matrix;
}

/* Function to transpose a matrix */
double **transposeMatrix(double **mat, int n, int d)
{
    int i, j;
    double **transposed;
    transposed = mallocMatrix(d, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            transposed[j][i] = mat[i][j];
        }
    }
    return transposed;
}

double **multiplyMatrices(double **mat1, int rows1, int cols1, double **mat2, int rows2, int cols2)
{
    int i;
    int j;
    int k;
    double **result;

    /* Check if the matrices can be multiplied */

    /* Allocate memory for the result matrix */
    result = mallocMatrix(rows1, cols2);

    /* Multiply the two matrices */
    for (i = 0; i < rows1; i++)
    {
        for (j = 0; j < cols2; j++)
        {
            result[i][j] = 0;
            for (k = 0; k < cols1; k++)
            {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

double calculationPerCell(double h, double wh, double hhth)
{
    double result;
    result = h * (1 - BETA + BETA * (wh / hhth));
    return result;
}

double **calculation(double **H, double **W_H, double **H_HT_H, int H_n, int H_d)
{
    int i;
    int j;
    double **new_H;
    new_H = mallocMatrix(H_n, H_d);
    for (i = 0; i < H_n; i++)
    {
        for (j = 0; j < H_d; j++)
        {
            new_H[i][j] = calculationPerCell(H[i][j], W_H[i][j], H_HT_H[i][j]);
        }
    }
    return new_H;
}

void copyMatrix(double **sourceMatrix, double **destinationMatrix, int rows, int cols)
{
    int i;
    int j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            destinationMatrix[i][j] = sourceMatrix[i][j];
        }
    }
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

    if (argc >= 3)
    {
        mode = argv[1];
        input_file_name = argv[2];
    }
    ArrayInfo datapointsstruct = read_file_to_array(input_file_name);
    datapoints = datapointsstruct.array;

    number_datapoints = datapointsstruct.rows;

    // if mode is equal to sym
    if (strcmp(mode, "sym") == 0)
    {
        outputmatrix = sym(datapoints, number_datapoints, dimensionOfVector);
    }
    else if (strcmp(mode, "ddg") == 0)
    {
        tempmatrix = sym(datapoints, number_datapoints, dimensionOfVector);
        outputmatrix = ddg(tempmatrix, number_datapoints);
        // free tempmatrix
        for (int i = 0; i < number_datapoints; i++)
        {
            free(tempmatrix[i]);
        }
    }
    else if (strcmp(mode, "norm") == 0)
    {
        tempmatrix = sym(datapoints, number_datapoints, dimensionOfVector);
        outputmatrix = norm(tempmatrix, number_datapoints);
        // free tempmatrix
        for (int i = 0; i < number_datapoints; i++)
        {
            free(tempmatrix[i]);
        }
    }
    // else if (strcmp(mode, "symnmf") == 0)
    // {

    //     K = atoi(argv[3]);
    //     printf("K = %d\n", K);
    //     tempmatrix =
    //         outputmatrix = symnmf(datapoints, K, number_datapoints);
    // }
    else
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    // free datapoints
    for (int i = 0; i < number_datapoints; i++)
    {
        free(datapoints[i]);
    }
    // print outputmatrix
    for (int i = 0; i < number_datapoints; i++)
    {
        for (int j = 0; j < number_datapoints; j++)
        {
            if (j != 0)
                printf(",");
            printf("%.4f", outputmatrix[i][j]);
        }
        printf("\n");
    }
    // free outputmatrix
    for (int i = 0; i < number_datapoints; i++)
    {
        free(outputmatrix[i]);
    }
}

/*endregion MAIN*/