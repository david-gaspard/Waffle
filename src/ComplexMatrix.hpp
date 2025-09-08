/****
 * @date Created on 2025-08-05 at 16:01:59 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the dense complex matrix class.
 ***/
#ifndef _COMPLEX_MATRIX_H
#define _COMPLEX_MATRIX_H
#include "RealMatrix.hpp"

/**
 * Class defining a dense complex matrix.
 */
class ComplexMatrix {
    
    private:
    
    dcomplex* data;  // Array containing the matrix elements in column-major order (LAPACK's convention).
    int nrow;       // Number of rows of the matrix.
    int ncol;       // Number of columns of the matrix.
    
    public:
    
    // Constructor/Destructor:
    ComplexMatrix(const int nrow, const int ncol);
    ComplexMatrix(const ComplexMatrix& a);
    ComplexMatrix(const RealMatrix& a);
    ~ComplexMatrix();
    
    // Getters/Setters:
    int getNrow() const;
    int getNcol() const;
    dcomplex& operator()(const int i, const int j) const;
    ComplexMatrix& operator=(const ComplexMatrix& a);
    
    // Mathematical operations:
    double norm() const;
    ComplexMatrix conj() const;
    RealMatrix real() const;
    RealMatrix imag() const;
    friend ComplexMatrix operator+(const ComplexMatrix& a, const ComplexMatrix& b);
    friend ComplexMatrix operator-(const ComplexMatrix& a, const ComplexMatrix& b);
    friend ComplexMatrix operator*(const dcomplex scalar, const ComplexMatrix& a);
    friend ComplexMatrix operator*(const ComplexMatrix& a, const dcomplex scalar);
    friend ComplexMatrix operator*(const ComplexMatrix& a, const ComplexMatrix& b);
    friend ComplexMatrix diagmul(const ComplexMatrix& a, const ComplexMatrix& b);
    friend void solve(const ComplexMatrix& a, const ComplexMatrix& b, ComplexMatrix& x);
    friend void svd(const ComplexMatrix& a, RealMatrix& s, ComplexMatrix& u, ComplexMatrix& v);
    friend ComplexMatrix apply(dcomplex(*f)(dcomplex), const ComplexMatrix& a);
    
    // Print methods:
    void print(const std::string& name) const;
    void savePNG(const std::string& filename) const;
    
};

// Constant matrices generators:
ComplexMatrix identityMatrix(const int n);
ComplexMatrix diagonalMatrix(const ComplexMatrix& diag, const int nrow, const int ncol);
ComplexMatrix gaussianRandomMatrix(const int nrow, const int ncol, const double sigma, const uint64_t seed);
double laplacianEigenvalue(const int i, const int n);
ComplexMatrix modalMatrix(const int n);
ComplexMatrix openingMatrix(const dcomplex kh2, const int n);
ComplexMatrix openingMatrix_old(const dcomplex kh2, const int n);

#endif
