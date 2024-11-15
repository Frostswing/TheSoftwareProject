#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

void printAndExitModule()
{
    printf("An Error Has Occurred");
    exit(1);
}

/* Python wrapper functions */
static PyObject *py_sym(PyObject *self, PyObject *args)
{
    PyObject *py_list;
    int n, m, i, j;
    double **X;
    double **result;
    PyObject *py_row;
    PyObject *py_val;
    PyObject *py_result;

    /* Parse Python arguments */
    if (!PyArg_ParseTuple(args, "Oii", &py_list, &n, &m))
    {
        return NULL;
    }

    /* Convert Python list to C array */
    X = malloc(n * sizeof(double *));
    if (X == NULL)
    {
        printAndExitModule();
    }
    for (i = 0; i < n; i++)
    {
        py_row = PyList_GetItem(py_list, i);
        X[i] = malloc(m * sizeof(double));
        if (X[i] == NULL)
        {
            printAndExitModule();
        }
        for (j = 0; j < m; j++)
        {
            py_val = PyList_GetItem(py_row, j);
            if (!py_val)
            {
                printAndExitModule();
            }
            X[i][j] = PyFloat_AsDouble(py_val);
        }
    }

    result = sym(X, n, m);

    /* Convert C array back to Python list */
    py_result = PyList_New(n);
    for (i = 0; i < n; i++)
    {
        py_row = PyList_New(n);
        for (j = 0; j < n; j++)
        {
            PyList_SetItem(py_row, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(py_result, i, py_row);
    }

    /* Free allocated C array and result array */
    for (i = 0; i < n; i++)
    {
        free(X[i]);
        free(result[i]);
    }
    free(X);
    free(result);

    return py_result;
}

static PyObject *py_ddg(PyObject *self, PyObject *args)
{
    PyObject *py_list;
    PyObject *py_row;
    PyObject *py_val;
    PyObject *py_result;
    int n, i, j;
    double **A;
    double **result;

    /* Parse Python arguments */
    if (!PyArg_ParseTuple(args, "Oi", &py_list, &n))
    {
        return NULL;
    }

    /* Convert Python list to C array */
    A = malloc(n * sizeof(double *));
    if (A == NULL)
    {
        printAndExitModule();
    }
    for (i = 0; i < n; i++)
    {
        py_row = PyList_GetItem(py_list, i);
        A[i] = malloc(n * sizeof(double));
        if (A[i] == NULL)
        {
            printAndExitModule();
        }
        for (j = 0; j < n; j++)
        {
            py_val = PyList_GetItem(py_row, j);
            A[i][j] = PyFloat_AsDouble(py_val);
        }
    }

    /* Call C function ddg */
    result = ddg(A, n);

    /* Convert C array to Python list */
    py_result = PyList_New(n);
    for (i = 0; i < n; i++)
    {
        py_row = PyList_New(n);
        for (j = 0; j < n; j++)
        {
            PyList_SetItem(py_row, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(py_result, i, py_row);
    }

    /* Free allocated C array and result array */
    for (i = 0; i < n; i++)
    {
        free(A[i]);
        free(result[i]);
    }
    free(A);
    free(result);

    return py_result;
}

static PyObject *py_norm(PyObject *self, PyObject *args)
{
    PyObject *py_list;
    PyObject *py_row;
    PyObject *py_val;
    int n, i, j;
    double **A;
    double **result;

    /* Parse Python arguments */
    if (!PyArg_ParseTuple(args, "Oi", &py_list, &n))
    {
        return NULL;
    }

    /* Convert Python list to C array */
    A = malloc(n * sizeof(double *));
    if (A == NULL)
    {
        printAndExitModule();
    }
    for (i = 0; i < n; i++)
    {
        py_row = PyList_GetItem(py_list, i);
        A[i] = malloc(n * sizeof(double));
        if (A[i] == NULL)
        {
            printAndExitModule();
        }
        for (j = 0; j < n; j++)
        {
            py_val = PyList_GetItem(py_row, j);
            A[i][j] = PyFloat_AsDouble(py_val);
        }
    }

    /* Call C function norm */
    result = norm(A, n);

    /* Convert C array to Python list */
    PyObject *py_result = PyList_New(n);
    for (i = 0; i < n; i++)
    {
        PyObject *py_row = PyList_New(n);
        for (j = 0; j < n; j++)
        {
            PyList_SetItem(py_row, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(py_result, i, py_row);
    }

    /* Free allocated C array and result array */
    for (i = 0; i < n; i++)
    {
        free(A[i]);
        free(result[i]);
    }
    free(A);
    free(result);

    return py_result;
}

static PyObject *py_symnmf(PyObject *self, PyObject *args)
{
    PyObject *py_list_H, *py_list_W;
    double **H, **W;
    int n, k, i, j;
    double **result;

    /* Parse Python arguments */
    if (!PyArg_ParseTuple(args, "OOii", &py_list_H, &py_list_W, &n, &k))
    {
        return NULL;
    }

    /* Convert Python list to C array for H */
    H = malloc(n * sizeof(double *));
    if (H == NULL)
    {
        printAndExitModule();
    }
    for (i = 0; i < n; i++)
    {
        PyObject *py_row = PyList_GetItem(py_list_H, i);
        H[i] = malloc(k * sizeof(double));
        if (H[i] == NULL)
        {
            printAndExitModule();
        }
        for (j = 0; j < k; j++)
        {
            PyObject *py_value = PyList_GetItem(py_row, j);
            H[i][j] = PyFloat_AsDouble(py_value);
        }
    }

    /* Convert Python list to C array for W */
    W = malloc(n * sizeof(double *));
    if(W == NULL)
    {
        printAndExitModule();
    }
    for (i = 0; i < n; i++)
    {
        PyObject *py_row = PyList_GetItem(py_list_W, i);
        W[i] = malloc(n * sizeof(double));
        if (W[i] == NULL)
        {
            printAndExitModule();
        }
        for (j = 0; j < n; j++)
        {
            PyObject *py_value = PyList_GetItem(py_row, j);
            W[i][j] = PyFloat_AsDouble(py_value);
        }
    }

    /* Call C function symnmf */
    result = symnmf(H, W, n, k);

    /* Convert C array to Python list */
    PyObject *py_result = PyList_New(n);
    for (i = 0; i < n; i++)
    {
        PyObject *py_row = PyList_New(k);
        for (j = 0; j < k; j++)
        {
            PyList_SetItem(py_row, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(py_result, i, py_row);
    }

    return py_result;
}

/* Define the methods that will be available in the Python module */
static PyMethodDef symnmfC_methods[] = {
    {"sym", py_sym, METH_VARARGS, "Calculate sym"},
    {"ddg", py_ddg, METH_VARARGS, "Calculate ddg"},
    {"norm", py_norm, METH_VARARGS, "Calculate norm"},
    {"symnmf", py_symnmf, METH_VARARGS, "Calculate symnmf"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfC_module = {
    PyModuleDef_HEAD_INIT,
    "symnmfC", /* name of module */
    NULL,      /* module documentation, may be NULL */
    -1,        /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables */
    symnmfC_methods
};

/* Initialization function for the module */
PyMODINIT_FUNC PyInit_symnmfC(void)
{
    return PyModule_Create(&symnmfC_module);
}
