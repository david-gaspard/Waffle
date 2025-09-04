/****
 * @date Created on 2025-08-29 at 13:16:10 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the SparseComplexMatrix class used to create and invert sparse complex matrices.
 * This version is based on a vector of triplets and a final sort to accelerate the matrix filling.
 ***/
#ifndef _SPARSE_COMPLEX_MATRIX_H
#define _SPARSE_COMPLEX_MATRIX_H
#include "ComplexMatrix.hpp"
#include <vector>
#include <iostream>

/**
 * Defines the element of a sparse complex matrix. This is known as the "triplet form".
 */
struct Triplet {
    
    int i;       // Row index of the element.
    int j;       // Column index of the element.
    dcomplex a;  // Matrix element.
    
};

std::ostream& operator<<(std::ostream& os, const Triplet& t);

/**
 * Defines the SparseComplexMatrix object.
 */
class SparseComplexMatrix {
    
    private: 
    
    std::vector<Triplet> triplet;  // Defines the list of triplets. Note that this vector is always sorted.
    int nrow;       // Number of rows of the sparse matrix.
    int ncol;       // Number of columns of the sparse matrix.
    bool symmetric; // Flag indicating if the matrix is symmetric or not (within a certain tolerance).
    bool sorted;    // Flag indicating if the sparse matrix is sorted by finalize() method.
    
    public:
    
    // Constructor/Destructor:
    SparseComplexMatrix(const int nrow, const int ncol);
    
    // Getters/Setters:
    int getNrow() const;
    int getNcol() const;
    int64_t getNnz() const;
    double density() const;
    void allocate(const int nnz);
    bool isSymmetric() const;
    dcomplex& operator()(const int i, const int j);
    dcomplex get(const int i, const int j) const;
    void checkSorted(const std::string& name) const;
    void finalize();
    
    // Print methods:
    void printSummary(const std::string& name) const;
    void print(const std::string& name) const;
    void saveImage(const std::string& filename) const;
    
    // Mathematical operations:
    double norm() const;
    SparseComplexMatrix conj() const;
    friend SparseComplexMatrix operator*(const dcomplex scalar, const SparseComplexMatrix& a);
    friend SparseComplexMatrix operator*(const SparseComplexMatrix& a, const dcomplex scalar);
    friend ComplexMatrix operator*(const SparseComplexMatrix& a, const ComplexMatrix& b);
    friend void solveUmfpack(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x);
    friend void solveMumps(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x);
    
    // Private computational methods:
    private:
    
    void checkIndices(const int i, const int j) const;
    void computeSymmetry();
    
};

#endif
