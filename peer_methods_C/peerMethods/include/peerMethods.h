/**
 * @author Vincenzo Iannucci
 * @brief The library provides an implementation for the main function
 * for solving peer method.
 * **/
#ifndef peerMethods_h
#define peerMethods_h

#ifdef __cplusplus
extern "C" {
#endif

/* Number of stages. */
#define STAGES 2

/* Variable that need to be set in the calling method. */
extern double a, B1, B2, F, H, S, d, D, L;
extern int M;

/*********************************************************************************
 * This struct has been created with the only purpose to return the value
 * obtained by the fPeerMethod function.
 ********************************************************************************/
typedef struct {
    double *yT;
    int yT_size;
    double *y;
    int y_rows;
    int y_cols;
    double *t;
    int t_size;
} return_values;

/**
 * @brief Initialize the struct return_values.
 * @param rv pointer to the struct return_values
*/
void initReturnStruct(return_values *rv);

/**
 * @brief Save the struct return_values in a file.
 * @param fileName the name of the file
 * @param rv pointer to the struct return
 * @return 0 if ok, 1 otherwise.
*/
int saveResultsInFile(const char* fileName, return_values result);

/**
 * @brief Build the matrix L.
 * This is an helping function that builds the matrix L.
 * @param L returning pointer to the matrix
 * @param LSize return the size of the matrix
 * @param Delta_x the value of the delta
*/
void computeLMatrix(double **L, int *LSize, double Delta_x);

/**
 * @brief Applies the Sherratt method.
 * @param y0 pointer to the y0 vector
 * @param y0Size size of the y0 vector
 * @param L pointer to the matrix L
 * @param LSize size of the matrix
 * @param sherrattSize returing size of the vector calculated by the function
 * @return a pointer the resulting vector after applying the Sherratt method.
*/
double *Sherratt(const double *y0, int y0Size, const double *L, int Lsize, int *sherrattSize);

/**
 * @brief Implicit fourth order method to solving ODE (Ordinary Differential Equation).
 * @param h number of conditions to achieve the solution
 * @param t0 starting time
 * @param y0 pointer to the y0 vector
 * @param y0Size size of the y0 vector
 * @param L pointer to the matrix L
 * @param LSize size of the matrix
 * @param ySize size of the result vector
 * @return a pointer to the y resulting vector.
*/
double *RungeKutta4th(double h, double t0, const double *y0, int y0Size, const double *L, int Lsize, int *ySize);

/**
 * @brief Compute the PDE (Partial Differential Equation) using the MOL (Method Of Lines).
 * The function computes the PDE (Partial Differential Equation) using MOL (Method Of Lines) and deriving 
 * a large system of ODE (Ordinary Differential Equation). Than, it solves the ODE system using the Runge
 * Kutta method of the fourth order.
 * @param N the size of the temporal grid
 * @param t_span an array representing the temporal grid itself
 * @param t_span_size the spatial dimension of the temporal grid
 * @param L pointer to the matrix L
 * @param LSize size of the matrix
 * @param y0 pointer to the y0 vector
 * @param y0Size size of the y0 vector
 * @param collect_result size of the result vector
 * @return a pointer to the y resulting vector.
*/
void fPeerClassic_twoStages(int N, double *t_span, int t_span_size, const double *L, int Lsize, const double *y0, int y0_size, return_values *collect_result);

/*************************************************
 *          Wrapper functions
 ***********************************************/

/**
 * @brief Function wrapper for malloc() function.
 * @param size Size of the memory allocated
 * @return a pointer to the allocated memory
*/
void *Malloc(size_t size);

/**
 * @brief Function wrapper for calloc() function.
 * @param nmemb number of elements to allocate
 * @param size Size of the memory allocated
 * @return a pointer to the allocated memory
*/
void *Calloc(size_t nmemb, size_t size);

/*************************************************
 *          Initialization functions
 ***********************************************/
void initializeRandomVector(double *vector, int N);
void initializeRandomMatrix(double *matrix, int M, int N);
int initMatrixByRowWithValuesFromVector(double *matrix, int M, int N, double *vector, int vector_size);
void initVectorWAnotherVector(double *newVector, double *oldVector, int n);

/*************************************************
 *  Free all the memory dynamically allocated
 ***********************************************/
void freeEverything(void *arg1, ...);

/*******************************************************
 * Return a vector in which the first element is 
 * first, the last element is last with a step that is 
 * represented by step parameters. N represent the
 * size of the array at the end.
 *****************************************************/
double *intervalDiscretization(double first, double last, double step, int *N);

/*******************************************************
 * Return a square matrix of size N in which the  
 * elemets at the position i == j assume value 1
 *****************************************************/
double *eyeD(int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 1.
 *****************************************************/
double *onesD(int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 0.
 *****************************************************/
double *zerosD(int N);

/*******************************************************
 * Return a matrix of size M * N in which every element
 * has value 0.
 *****************************************************/
double *zerosMatrixD(int M, int N);

/*******************************************************
 * Return a matrix of matrix_size x matrix_size elements 
 * in which each element of the k-th diagonal is an 
 * element of the vector passed by parameter.
 *****************************************************/
double *diagD(double *vector, int size, int k, int *matrix_size);

/*******************************************************
 * Packing three matrices side by side into one. It's
 * an utility function for the threeBlockDiagD function
 *****************************************************/
double *packThreeMatrices(int n, double *A, double *B, double *C);

/*******************************************************
 * Return a block diag matrix of three matrices passed
 * by arguments of size n x n.
 *****************************************************/
double *threeBlockDiagD(int n, double *A, double *B, double *C, int *blckSize);

/*******************************************************
 * Packing three vectors side by side into one.
 *****************************************************/
double *packThreeVectors(int n, double *A, double *B, double *C, int *newDimension);

/****************************************************
 * Generates n points. The spacing between the points 
 * is (x2-x1)/(n-1)
 ***************************************************/
double *linspace(double x1, double x2, int n);

#ifdef __cplusplus
}
#endif
#endif // !peerMethods_h