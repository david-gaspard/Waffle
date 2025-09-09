/****
 * @date Created on 2025-08-12 at 15:54:25 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the dense real matrix class.
 ***/
#ifndef _REAL_MATRIX_H
#define _REAL_MATRIX_H
#include "Constants.hpp"

class ComplexMatrix;

/**
 * Class defining a dense real matrix.
 */
class RealMatrix {
    
    friend ComplexMatrix;  // State that ComplexMatrix may access the private fields of RealMatrix.
    
    private:
    
    double* data;  // Array containing the matrix elements in column-major order (LAPACK's convention).
    int nrow;    // Number of rows of the matrix.
    int ncol;    // Number of columns of the matrix.
    
    public:
    
    // Constructor/Destructor:
    RealMatrix(const int nrow, const int ncol);
    RealMatrix(const RealMatrix& a);
    ~RealMatrix();
    
    // Getters/Setters:
    int getNrow() const;
    int getNcol() const;
    double& operator()(const int i, const int j) const;
    RealMatrix& operator=(const RealMatrix& a);
    
    // Mathematical operations:
    double norm() const;
    double sum() const;
    double max() const;
    double min() const;
    RealMatrix transpose() const;
    friend RealMatrix operator+(const RealMatrix& a, const RealMatrix& b);
    friend RealMatrix operator-(const RealMatrix& a, const RealMatrix& b);
    friend RealMatrix& operator+=(RealMatrix& a, const RealMatrix& b);
    friend RealMatrix operator*(const double scalar, const RealMatrix& a);
    friend RealMatrix operator*(const RealMatrix& a, const double scalar);
    friend RealMatrix operator*(const RealMatrix& a, const RealMatrix& b);
    friend void svd(const ComplexMatrix& a, RealMatrix& s, ComplexMatrix& u, ComplexMatrix& v);
    
    // Print methods:
    void print(const std::string& name) const;
    
};

#endif
