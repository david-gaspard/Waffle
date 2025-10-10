/****
 * @date Created on 2025-08-12 at 15:58:10 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the methods for dense real matrices.
 ***/
#include "RealMatrix.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>

/**
 * Constructor of the dense matrix object.
 */
RealMatrix::RealMatrix(const int nrow, const int ncol) {
    if (nrow <= 0) {
        throw std::invalid_argument("In RealMatrix(): Number of rows cannot be negative or zero.");
    }
    if (ncol <= 0) {
        throw std::invalid_argument("In RealMatrix(): Number of columns cannot be negative or zero.");
    }
    this->nrow = nrow;
    this->ncol = ncol;
    data = new double[nrow*ncol]();  // Allocate space for the matrix (column-major format) with zero initialization.
}

/**
 * Explicit copy constructor because of the underlying pointer "data".
 */
RealMatrix::RealMatrix(const RealMatrix& a) {
    nrow = a.nrow;
    ncol = a.ncol;
    const int nt = nrow*ncol;
    data = new double[nt];
    std::copy(a.data, a.data + nt, data);
}

/**
 * Destructor of a matrix. Frees the memory allocated by the constructor.
 */
RealMatrix::~RealMatrix() {
    delete[] data;
}


/***************************************************************************************************
 * GETTERS AND SETTERS
 ***************************************************************************************************/

/**
 * Returns the number of rows of this matrix.
 */
int RealMatrix::getNrow() const {
    return nrow;
}

/**
 * Returns the number of columns of this matrix.
 */
int RealMatrix::getNcol() const {
    return ncol;
}

/**
 * Returns a reference to the element (i, j) of the matrix.
 * The row index is "i", and the column index is "j".
 * Note that the indices "i" and "j" are based to zero (C/C++ convention).
 */
double& RealMatrix::operator()(const int i, const int j) const {
    if (i >= nrow || j >= ncol || i < 0 || j < 0) {
        std::string msg = "In Matrix(): Invalid index (i=" + std::to_string(i) + ", j=" + std::to_string(j) 
                        + ") given nrow=" + std::to_string(nrow) + " and ncol=" + std::to_string(ncol) + ".";
        throw std::out_of_range(msg);
    }
    return data[i + j*nrow];
}

/**
 * Overloads the assignement operator to perform a deep copy.
 */
RealMatrix& RealMatrix::operator=(const RealMatrix& a) {
    if (this == &a) {
        throw std::invalid_argument("In RealMatrix operator=(): Invalid self assignment.");
    }
    if (nrow != a.nrow || ncol != a.ncol) {
        std::string msg = "In RealMatrix operator=(): Matrices have different sizes. LHS = (" 
                         + std::to_string(nrow) + ", " + std::to_string(ncol) 
                         + "), RHS = (" + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    std::copy(a.data, a.data + nrow*ncol, data); // Deep copy of the array (should be faster than loop).
    
    return *this;
}

/***************************************************************************************************
 * ALGEBRAIC OPERATIONS ON MATRICES
 ***************************************************************************************************/

/**
 * Computes the Frobenius norm of the present matrix.
 * This is equivalent to the 2-norm for vectors.
 */
double RealMatrix::norm() const {
    
    double a, sumsq = 0.;
    
    for (int l = 0; l < nrow*ncol; l++) {// Loop over all elements.
        a = data[l];
        sumsq += a*a;
    }
    return std::sqrt(sumsq);
}

/**
 * Computes the sum of all matrix elements of the present matrix.
 */
double RealMatrix::sum() const {
    
    //double sum = 0.;
    //
    //for (int l = 0; l < nrow*ncol; l++) {// Loop over all elements.
    //    sum += data[l];
    //}
    //return sum;
    
    return std::reduce(data, data + nrow*ncol);
}

/**
 * Computes the mean of all matrix elements of the present matrix.
 */
double RealMatrix::mean() const {
    return sum()/(static_cast<double>(nrow) * static_cast<double>(ncol));
}

/**
 * Computes the standard deviation of all matrix elements of the present matrix.
 */
double RealMatrix::stddev() const {
    
    const double avg = mean();
    double delta, sumsq = 0.;
    
    for (int l = 0; l < nrow*ncol; l++) {// Loop over all elements.
        delta = data[l] - avg;
        sumsq += delta * delta;
    }
    
    return std::sqrt(sumsq/(static_cast<double>(nrow) * static_cast<double>(ncol)));
}

/**
 * Returns the maximum value of all matrix elements of the present matrix.
 */
double RealMatrix::max() const {
    return *std::max_element(data, data+nrow*ncol);
}

/**
 * Returns the maximum value of all matrix elements of the present matrix.
 */
double RealMatrix::min() const {
    return *std::min_element(data, data+nrow*ncol);
}

/**
 * Returns the transpose of the present matrix.
 */
RealMatrix RealMatrix::transpose() const {
    
    RealMatrix c(ncol, nrow);
    
    for (int j = 0; j < ncol; j++) {// Loop over the columns.
        for (int i = 0; i < nrow; i++) {// Loop over the rows.
            c.data[j + i*ncol] = data[i + j*nrow];
        }
    }
    
    return c;
}

/**
 * Defines the addition of two matrices of same dimensions.
 */
RealMatrix operator+(const RealMatrix& a, const RealMatrix& b) {
    
    if (a.nrow != b.nrow || a.ncol != b.ncol) {// Check for possible incompatibilities.
        std::string msg = "In operator+(): Incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    RealMatrix sum(a.nrow, a.ncol);
    
    for (int l = 0; l < a.nrow*a.ncol; l++) {
        sum.data[l] = a.data[l] + b.data[l];
    }
    
    return sum;
}

/**
 * Defines the subtraction of two matrices of same dimensions.
 */
RealMatrix operator-(const RealMatrix& a, const RealMatrix& b) {
    
    if (a.nrow != b.nrow || a.ncol != b.ncol) {// Check for possible incompatibilities.
        std::string msg = "In operator-(): Incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    RealMatrix diff(a.nrow, a.ncol);
    
    for (int l = 0; l < a.nrow*a.ncol; l++) {
        diff.data[l] = a.data[l] - b.data[l];
    }
    
    return diff;
}

/**
 * Overload the increment operator.
 */
RealMatrix& operator+=(RealMatrix& a, const RealMatrix& b) {
    
    if (a.nrow != b.nrow || a.ncol != b.ncol) {// Check for possible incompatibilities.
        std::string msg = "In operator+=(): Incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    for (int l = 0; l < a.nrow*a.ncol; l++) {
        a.data[l] += b.data[l];
    }
    
    return a;
}

/**
 * Multiplication of a matrix by a scalar.
 */
RealMatrix operator*(const double scalar, const RealMatrix& a) {
    
    RealMatrix b(a.nrow, a.ncol);
    
    for (int l = 0; l < a.nrow*a.ncol; l++) {// Loop on the matrix elements.
        b.data[l] = scalar * a.data[l];
    }
    
    return b;
}

/**
 * Multiplication of a matrix by a scalar the other way around.
 */
RealMatrix operator*(const RealMatrix& a, const double scalar) {
    return scalar*a;
}

/**
 * Product of double-precision real matrices using BLAS.
 */
extern "C" void dgemm_(const char& transa, const char& transb, const int& m, const int& n, const int& k, const double& alpha, 
                       const double* a, const int& lda, const double* b, const int& ldb, const double& beta, double* c, const int& ldc);

/**
 * Defines the product of two matrices of compatible dimensions.
 * This function is a user-friendly wrapper to dgemm_().
 */
RealMatrix operator*(const RealMatrix& a, const RealMatrix& b) {
    
    if (a.ncol != b.nrow) {// Check for possible incompatibilities.
        std::string msg = "In operator*(): Product of incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    RealMatrix c(a.nrow, b.ncol);
    
    const char trans = 'N';   // Indicate that the matrices must not be transposed.
    const double alpha = 1.;  // Constant "alpha" in C = alpha*A*B + beta*C (see "dgemm" in BLAS manual).
    const double beta = 0.;   // Constant "beta" in C = alpha*A*B + beta*C (see "dgemm" in BLAS manual).
    
    dgemm_(trans, trans, a.nrow, b.ncol, a.ncol, alpha, a.data, a.nrow, b.data, b.nrow, beta, c.data, a.nrow);
    
    return c;
}

/***************************************************************************************************
 * PRINTING METHODS
 ***************************************************************************************************/

/**
 * Print the matrix to standard output.
 */
void RealMatrix::print(const std::string& name) const {
    
    std::cout << TAG_INFO << name << " = [\n";
    
    for (int i = 0; i < nrow; i++) {// Loop over the rows.
        std::cout << "\t[";
        for (int j = 0; j < ncol; j++) {// Loop over the columns.
            std::cout << " " << data[i + j*nrow] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "]" << std::endl;
    
}
