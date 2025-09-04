/****
 * @date Created on 2025-08-04 at 17:38:31 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the SparseComplexMatrix class.
 ***/
#include "SparseComplexMatrix.hpp"
#include <iostream>

/**
 * Print generic informations about the matrix.
 */
void printInfoMatrix(const std::string& name, const SparseComplexMatrix& a) {
    a.printSummary(name);
    a.print(name);
}

/**
 * Test the filling of the SparseComplexMatrix.
 */
int testFilling1() {
    std::cout << "====== SPARSE COMPLEX MATRIX - TEST FILLING #1 ======" << std::endl;
    
    SparseComplexMatrix matrix(5, 3);
    
    matrix(2, 1) = dcomplex(1., 2.);
    matrix(4, 0) = dcomplex(3., 4.);
    matrix.finalize();
    
    std::cout << TAG_INFO << "matrix.getNrow = " << matrix.getNrow() << "\n";
    std::cout << TAG_INFO << "matrix.getNcol = " << matrix.getNcol() << "\n";
    std::cout << TAG_INFO << "matrix.getNnz  = " << matrix.getNnz() << "\n";
    std::cout << TAG_INFO << "matrix(2, 2)   = " << matrix.get(2, 2) << "\n";
    std::cout << TAG_INFO << "matrix(2, 1)   = " << matrix.get(2, 1) << "\n";
    std::cout << TAG_INFO << "matrix(4, 0)   = " << matrix.get(4, 0) << "\n";
    std::cout << TAG_INFO << "matrix.getNnz  = " << matrix.getNnz() << "\n";
    
    printInfoMatrix("Matrix", matrix);
    
    //printInfoMatrix("Matrix*", matrix.conj());
    
    return 0;
}

/**
 * Test the filling of the SparseComplexMatrix.
 */
int testFilling2() {
    std::cout << "====== SPARSE COMPLEX MATRIX - TEST FILLING #2 ======" << std::endl;
    
    const int n = 5;
    dcomplex elem;
    SparseComplexMatrix matrix(n, n);
    
    for (int i = 0; i < n; i++) {
        matrix(i, i) = dcomplex(3*i, -i);
        for (int j = 0; j < i; j++) {// Loop over the lower triangle.
            dcomplex elem = dcomplex(i + 2*j, j - 2*i);
            matrix(i, j) = elem;
            matrix(j, i) = elem;
        }
    }
    matrix.finalize();
    printInfoMatrix("Matrix", matrix);
    
    return 0;
}

/**
 * Test the hybrid product of a sparse matrix by a dense matrix.
 */
int testHybridProduct() {
    std::cout << "====== TEST HYBRID PRODUCT SPARSE*DENSE ======" << std::endl;
    
    SparseComplexMatrix a(5, 3);
    a(0, 0) = 2;   a(1, 0) = 3;  a(0, 1) = 3;
    a(2, 1) = -1;  a(4, 1) = 4;  a(1, 2) = 4;
    a(2, 2) = -3;  a(3, 2) = 1;  a(4, 2) = 2;
    a.finalize();
    
    ComplexMatrix b(3, 2);
    b(0, 0) = 1;  b(0, 1) = 4;
    b(1, 0) = 2;  b(1, 1) = 5;
    b(2, 0) = 3;  b(2, 1) = 6;
    
    ComplexMatrix c_expc(5, 2);
    c_expc(0, 0) =   8;  c_expc(0, 1) =  23;
    c_expc(1, 0) =  15;  c_expc(1, 1) =  36;
    c_expc(2, 0) = -11;  c_expc(2, 1) = -23;
    c_expc(3, 0) =   3;  c_expc(3, 1) =   6;
    c_expc(4, 0) =  14;  c_expc(4, 1) =  32;
    
    ComplexMatrix c = a*b;
    
    //c.print("C = A*B");
    
    std::cout << TAG_INFO << "Error = " << (c - c_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the solution of a sparse linear system A.x = b using the UMFPACK solver.
 */
int testSolveUmfpack1() {
    std::cout << "====== TEST SOLVE UMFPACK #1 - NONSYMMETRIC MATRIX ======" << std::endl;
    
    const int n = 5;
    SparseComplexMatrix a(n, n);
    a(0, 0) = 2;   a(1, 0) = 3;  a(0, 1) = 3;  a(2, 1) = -1;  a(4, 1) = 4;  a(1, 2) = 4;
    a(2, 2) = -3;  a(3, 2) = 1;  a(4, 2) = 2;  a(2, 3) = 2;   a(1, 4) = 6;  a(4, 4) = 1;
    a.finalize();
    
    //printInfoMatrix("A", a);
    
    SparseComplexMatrix b(n, 1);
    b(0, 0) = 8;  b(1, 0) = 45;  b(2, 0) = -3;  b(3, 0) = 3;  b(4, 0) = 19; 
    b.finalize();
    
    ComplexMatrix x(n, 1), x_expc(n, 1);
    x_expc(0, 0) = 1;  x_expc(1, 0) = 2;  x_expc(2, 0) = 3;  x_expc(3, 0) = 4;  x_expc(4, 0) = 5;
    
    solveUmfpack(a, b, x);
    
    //x.print("Solution x");
    
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the solution of a sparse linear system A.x = b using the UMFPACK solver.
 * This tests the case of a symmetric matrix.
 */
int testSolveUmfpack2() {
    std::cout << "====== TEST SOLVE UMFPACK #2 - SYMMETRIC MATRIX ======" << std::endl;
    
    const int n = 5;
    SparseComplexMatrix a(n, n);
    a(0, 0) = 2;  a(1, 0) = 6;  a(0, 1) = 6;  a(2, 1) =  3;  a(4, 1) = 10;  a(1, 2) = 3;  a(2, 2) = -3;
    a(3, 2) = 3;  a(4, 2) = 2;  a(2, 3) = 3;  a(1, 4) = 10;  a(2, 4) =  2;  a(4, 4) = 1;
    a.finalize();
    
    //printInfoMatrix("A", a);
    
    SparseComplexMatrix b(n, 1);
    b(0, 0) = 14;  b(1, 0) = 65;  b(2, 0) = 19;  b(3, 0) = 9;  b(4, 0) = 31;
    b.finalize();
    
    ComplexMatrix x(n, 1), x_expc(n, 1);
    x_expc(0, 0) = 1;  x_expc(1, 0) = 2;  x_expc(2, 0) = 3;  x_expc(3, 0) = 4;  x_expc(4, 0) = 5;
    
    solveUmfpack(a, b, x);
    
    //x.print("Solution x");
    
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the solution of a sparse linear system A.xb = b using the UMFPACK solver.
 * Case fully complex with multiple right-hand sides.
 */
int testSolveUmfpack3() {
    std::cout << "====== TEST SOLVE UMFPACK #3 - MULTIPLE RIGHT-HAND SIDES ======" << std::endl;
    
    // 1. Prepare matrix 'A':
    const int n = 5;
    SparseComplexMatrix a(n, n);
    a(0, 0) = dcomplex(2, -1);
    a(0, 1) = dcomplex(3, 2);
    a(1, 2) = dcomplex(5, -7);
    a(1, 4) = dcomplex(3, 2);
    a(2, 1) = dcomplex(-1, -1);
    a(2, 2) = dcomplex(-3, -2);
    a(2, 3) = dcomplex(2, 1);
    a(3, 0) = dcomplex(0, 1);
    a(3, 2) = dcomplex(1, 2);
    a(4, 1) = dcomplex(4, 3);
    a(4, 2) = dcomplex(-2, -1);
    a(4, 4) = dcomplex(-1, 5);
    a.finalize();
    
    //printInfoMatrix("A", a);
    
    // 2. Prepare independent term 'b':
    SparseComplexMatrix b(n, 3);
    b(0, 0) = dcomplex(  8,   3);  b(0, 1) = dcomplex( 33,   8);  b(0, 2) = dcomplex(  2,  -1);
    b(1, 0) = dcomplex( 30, -11);  b(1, 1) = dcomplex( 70, -36);  b(1, 2) = dcomplex(  8,  -5);
    b(2, 0) = dcomplex( -3,  -4);  b(2, 1) = dcomplex(-13, -14);  b(2, 2) = dcomplex( -3,  -2);
    b(3, 0) = dcomplex(  3,   7);  b(3, 1) = dcomplex(  8,  22);  b(3, 2) = dcomplex(  1,   3);
    b(4, 0) = dcomplex( -3,  28);  b(4, 1) = dcomplex(  2,  63);  b(4, 2) = dcomplex( -3,   4);
    b.finalize();
    
    // 3. Prepare expected solution 'x_expc':
    ComplexMatrix x(n, 3), x_expc(n, 3);
    x_expc(0, 0) = 1;  x_expc(0, 1) =  6;  x_expc(0, 2) = 1;
    x_expc(1, 0) = 2;  x_expc(1, 1) =  7;  x_expc(1, 2) = 0;
    x_expc(2, 0) = 3;  x_expc(2, 1) =  8;  x_expc(2, 2) = 1;
    x_expc(3, 0) = 4;  x_expc(3, 1) =  9;  x_expc(3, 2) = 0;
    x_expc(4, 0) = 5;  x_expc(4, 1) = 10;  x_expc(4, 2) = 1;
    
    solveUmfpack(a, b, x);
    
    //x.print("Solution x");
    
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the solution of a sparse linear system A.x = b using the MUMPS solver.
 */
int testSolveMumps1() {
    std::cout << "====== TEST SOLVE MUMPS #1 - NONSYMMETRIC MATRIX ======" << std::endl;
    
    const int n = 5;
    SparseComplexMatrix a(n, n);
    a(0, 0) = 2;   a(1, 0) = 3;  a(0, 1) = 3;  a(2, 1) = -1;  a(4, 1) = 4;  a(1, 2) = 4;
    a(2, 2) = -3;  a(3, 2) = 1;  a(4, 2) = 2;  a(2, 3) = 2;   a(1, 4) = 6;  a(4, 4) = 1;
    a.finalize();
    
    //printInfoMatrix("A", a);
    
    SparseComplexMatrix b(n, 1);
    b(0, 0) = 8;  b(1, 0) = 45;  b(2, 0) = -3;  b(3, 0) = 3;  b(4, 0) = 19; 
    b.finalize();
    
    ComplexMatrix x(n, 1), x_expc(n, 1);
    x_expc(0, 0) = 1;  x_expc(1, 0) = 2;  x_expc(2, 0) = 3;  x_expc(3, 0) = 4;  x_expc(4, 0) = 5;
    
    solveMumps(a, b, x);
    
    //x.print("Solution x");
    
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the solution of a sparse linear system A.x = b using the MUMPS solver.
 * This tests the case of a symmetric matrix.
 */
int testSolveMumps2() {
    std::cout << "====== TEST SOLVE MUMPS #2 - SYMMETRIC MATRIX ======" << std::endl;
    
    const int n = 5;
    SparseComplexMatrix a(n, n);
    a(0, 0) = 2;  a(1, 0) = 6;  a(0, 1) = 6;  a(2, 1) =  3;  a(4, 1) = 10;  a(1, 2) = 3;  a(2, 2) = -3;
    a(3, 2) = 3;  a(4, 2) = 2;  a(2, 3) = 3;  a(1, 4) = 10;  a(2, 4) =  2;  a(4, 4) = 1;
    a.finalize();
    
    //printInfoMatrix("A", a);
    
    SparseComplexMatrix b(n, 1);
    b(0, 0) = 14;  b(1, 0) = 65;  b(2, 0) = 19;  b(3, 0) = 9;  b(4, 0) = 31;
    b.finalize();
    
    ComplexMatrix x(n, 1), x_expc(n, 1);
    x_expc(0, 0) = 1;  x_expc(1, 0) = 2;  x_expc(2, 0) = 3;  x_expc(3, 0) = 4;  x_expc(4, 0) = 5;
    
    solveMumps(a, b, x);
    
    //x.print("Solution x");
    
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the solution of a sparse linear system A.xb = b using the MUMPS solver.
 * Case fully complex with multiple right-hand sides.
 */
int testSolveMumps3() {
    std::cout << "====== TEST SOLVE MUMPS #3 - MULTIPLE RIGHT-HAND SIDES ======" << std::endl;
    
    // 1. Prepare matrix 'A':
    const int n = 5;
    SparseComplexMatrix a(n, n);
    a(0, 0) = dcomplex(2, -1);
    a(0, 1) = dcomplex(3, 2);
    a(1, 2) = dcomplex(5, -7);
    a(1, 4) = dcomplex(3, 2);
    a(2, 1) = dcomplex(-1, -1);
    a(2, 2) = dcomplex(-3, -2);
    a(2, 3) = dcomplex(2, 1);
    a(3, 0) = dcomplex(0, 1);
    a(3, 2) = dcomplex(1, 2);
    a(4, 1) = dcomplex(4, 3);
    a(4, 2) = dcomplex(-2, -1);
    a(4, 4) = dcomplex(-1, 5);
    a.finalize();
    
    //printInfoMatrix("A", a);
    
    // 2. Prepare independent term 'b':
    SparseComplexMatrix b(n, 3);
    b(0, 0) = dcomplex(  8,   3);  b(0, 1) = dcomplex( 33,   8);  b(0, 2) = dcomplex(  2,  -1);
    b(1, 0) = dcomplex( 30, -11);  b(1, 1) = dcomplex( 70, -36);  b(1, 2) = dcomplex(  8,  -5);
    b(2, 0) = dcomplex( -3,  -4);  b(2, 1) = dcomplex(-13, -14);  b(2, 2) = dcomplex( -3,  -2);
    b(3, 0) = dcomplex(  3,   7);  b(3, 1) = dcomplex(  8,  22);  b(3, 2) = dcomplex(  1,   3);
    b(4, 0) = dcomplex( -3,  28);  b(4, 1) = dcomplex(  2,  63);  b(4, 2) = dcomplex( -3,   4);
    b.finalize();
    
    // 3. Prepare expected solution 'x_expc':
    ComplexMatrix x(n, 3), x_expc(n, 3);
    x_expc(0, 0) = 1;  x_expc(0, 1) =  6;  x_expc(0, 2) = 1;
    x_expc(1, 0) = 2;  x_expc(1, 1) =  7;  x_expc(1, 2) = 0;
    x_expc(2, 0) = 3;  x_expc(2, 1) =  8;  x_expc(2, 2) = 1;
    x_expc(3, 0) = 4;  x_expc(3, 1) =  9;  x_expc(3, 2) = 0;
    x_expc(4, 0) = 5;  x_expc(4, 1) = 10;  x_expc(4, 2) = 1;
    
    solveMumps(a, b, x);
    
    //x.print("Solution x");
    
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm() << "\n";
    
    return 0;
}

/**
 * Main function of the test of the SparseComplexMatrix class.
 */
int main(int argc, char** argv) {
    
    //testFilling1();
    //testFilling2();
    //testHybridProduct();
    testSolveUmfpack1();
    testSolveUmfpack2();
    testSolveUmfpack3();
    testSolveMumps1();
    testSolveMumps2();
    testSolveMumps3();
    
    return 0;
}
