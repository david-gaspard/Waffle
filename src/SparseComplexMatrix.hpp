/****
 * @date Created on 2025-08-27 at 14:43:49 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the SparseComplexMatrix class used to create and invert sparse complex matrices.
 * This object serves as an interface to UMFPACK and possibly other sparse linear solvers.
 ***/
#ifndef _SPARSE_COMPLEX_MATRIX_H
#define _SPARSE_COMPLEX_MATRIX_H
#include "ComplexMatrix.hpp"
#include <map>

/**
 * Defines the indices of the sparse matrix.
 * Note that the indices (i, j) must be between (0, 0) and (nrow-1, ncol-1).
 */
struct Indices {
    
    int i;   // Row index of the element.
    int j;   // Column index of the element.
    
};

/**
 * Defines the SparseComplexMatrix object.
 */
class SparseComplexMatrix {
    
    private: 
    
    std::map<Indices, dcomplex> data;  // List of matrix elements stored in column-major ordering.
        // Note: Maps ensure that search, insertion, and removal operations have logarithmic complexity O(log(N)),
        // and are thus much more efficient than Vectors for very large N. Maps greatly improve the construction speed of the sparse matrix.
    int nrow;       // Number of rows of the sparse matrix.
    int ncol;       // Number of columns of the sparse matrix.
    
    public:
    
    // Constructor/Destructor:
    SparseComplexMatrix(const int nrow, const int ncol);
    
    // Getters/Setters:
    int getNrow() const;
    int getNcol() const;
    int64_t getNnz() const;
    dcomplex& operator()(const int i, const int j);
    dcomplex get(const int i, const int j) const;
    double density() const;
    
    // Mathematical operations:
    double norm() const;
    SparseComplexMatrix conj() const;
    bool isSymmetric() const;
    friend SparseComplexMatrix operator*(const dcomplex scalar, const SparseComplexMatrix& a);
    friend SparseComplexMatrix operator*(const SparseComplexMatrix& a, const dcomplex scalar);
    friend ComplexMatrix operator*(const SparseComplexMatrix& a, const ComplexMatrix& b);
    friend void solveUmfpack(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x);
    
    // Print methods:
    void summary(const std::string& name) const;
    void print(const std::string& name) const;
    void saveImage(const std::string& filename) const;
    
    // Private computational methods:
    private:
    
    void checkIndices(const int i, const int j) const;
    
};

#endif
