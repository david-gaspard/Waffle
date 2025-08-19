/****
 * @date Created on 2025-08-04 at 13:43:20 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the SparseComplexMatrix class used to create and invert sparse coplex matrices.
 * This object serves as an interface to UMFPACK and MUMPS.
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
    
    int i;   // Row index of the element.
    int j;   // Column index of the element.
    dcomplex a;  // Matrix element.
    
};

bool compareTriplet(const Triplet& t1, const Triplet& t2);
std::ostream& operator<<(std::ostream& os, const Triplet& t);

/**
 * Defines the SparseComplexMatrix object.
 */
class SparseComplexMatrix {
    
    private: 
    
    std::vector<Triplet> triplet;  // Defines the list of triplets. Note that this vector is always sorted.
    int nrow;       // Number of rows of the sparse matrix.
    int ncol;       // Number of columns of the sparse matrix.
    int symmetry;   // Flag indicating if the matrix is symmetric (invariant under matrix transpose). 0=Not symmetric, 1=Symmetric.
                    // If this flag is 1, only the upper or lower part of the matrix is needed.
    
    public:
    
    // Constructor/Destructor:
    SparseComplexMatrix(const int nrow, const int ncol);
    
    // Getters:
    int getNrow() const;
    int getNcol() const;
    int64_t getNnz() const;
    int getSymmetry() const;
    dcomplex get(const int i, const int j) const;
    double density() const;
    
    // Setters:
    void setSymmetry(const int symmetry);
    dcomplex& operator()(const int i, const int j);
    
    // Mathematical operations:
    double norm() const;
    SparseComplexMatrix conj() const;
    bool isSymmetric() const;
    friend SparseComplexMatrix operator*(const dcomplex scalar, const SparseComplexMatrix& a);
    friend SparseComplexMatrix operator*(const SparseComplexMatrix& a, const dcomplex scalar);
    friend ComplexMatrix operator*(const SparseComplexMatrix& a, const ComplexMatrix& b);
    friend void solveUmfpack(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x);
    
    // Print methods:
    void printInfo(const std::string& name) const;
    void print(const std::string& name) const;
    void printImage(const std::string& filename) const;
    
    // Private computational methods:
    private:
    
    void checkIndices(const int i, const int j) const;
    friend void solveUmfpack_old(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x);
    friend void solveMumps(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x);
    
};

#endif
