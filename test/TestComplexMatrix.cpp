/****
 * @date Created on 2025-08-05 at 21:36:35 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the ComplexMatrix class.
 ***/
#include "ComplexMatrix.hpp"
#include <chrono>
#include <iostream>

/**
 * Test the Hermitian conjugation.
 */
int testHermitianConjugate() {
    std::cout << "====== TEST HERMITIAN CONJUGATION ======" << std::endl;
    
    ComplexMatrix a(4, 3);
    
    a(0, 0) = dcomplex(1,  2);  a(0, 1) = dcomplex( 3, -4);  a(0, 2) = dcomplex(-5, 6);  
    a(1, 0) = dcomplex(3, -1);  a(1, 1) = dcomplex(-3,  1);  a(1, 2) = dcomplex( 1, 1);  
    a(2, 0) = dcomplex(2,  0);  a(2, 1) = dcomplex(-1, -1);  a(2, 2) = dcomplex( 0, 1);  
    a(3, 0) = dcomplex(0, -2);  a(3, 1) = dcomplex( 4,  1);  a(3, 2) = dcomplex(10, 3);  
    
    a.print("Matrix A");
    
    a.conj().print("Matrix A*");
    
    return 0;
}

/**
 * Test the product of two rectangular complex matrices using BLAS's zgemm() routine.
 */
int testMatrixProduct() {
    std::cout << "====== TEST MATRIX PRODUCT ======" << std::endl;
    
    ComplexMatrix a(4, 3), b(3, 2), c_expc(4, 2);
    
    a(0, 0) = dcomplex(1,  2);  a(0, 1) = dcomplex( 3, -4);  a(0, 2) = dcomplex(-5, 6);  
    a(1, 0) = dcomplex(3, -1);  a(1, 1) = dcomplex(-3,  1);  a(1, 2) = dcomplex( 1, 1);  
    a(2, 0) = dcomplex(2,  0);  a(2, 1) = dcomplex(-1, -1);  a(2, 2) = dcomplex( 0, 1);  
    a(3, 0) = dcomplex(0, -2);  a(3, 1) = dcomplex( 4,  1);  a(3, 2) = dcomplex(10, 3);  
    
    b(0, 0) = dcomplex(5, -2);  b(0, 1) = dcomplex( 3, 4);
    b(1, 0) = dcomplex(0, -1);  b(1, 1) = dcomplex( 2, 1);
    b(2, 0) = dcomplex(1,  6);  b(2, 1) = dcomplex(-7, 5);
    
    c_expc(0, 0) = dcomplex(-36, -19); c_expc(0, 1) = dcomplex( 10, -62);
    c_expc(1, 0) = dcomplex(  9,  -1); c_expc(1, 1) = dcomplex( -6,   6);
    c_expc(2, 0) = dcomplex(  3,  -2); c_expc(2, 1) = dcomplex(  0,  -2);
    c_expc(3, 0) = dcomplex(-11,  49); c_expc(3, 1) = dcomplex(-70,  29);
    
    std::cout << TAG_INFO << "Error (A*B - C_expc) = " << (a * b - c_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the multiplication by a diagonal matrix.
 */
int testDiagmul() {
    std::cout << "====== TEST DIAGMUL ======" << std::endl;
    
    const int n1 = 5;
    const int n2 = 3;
    const double sigma = 1.;
    uint64_t seed = 1;
    
    ComplexMatrix a = gaussianRandomMatrix(n1, n2, sigma, seed);
    ComplexMatrix dl = gaussianRandomMatrix(n1, 1, sigma, seed+1);
    ComplexMatrix dr = gaussianRandomMatrix(1, n2, sigma, seed+2);
    
    ComplexMatrix cl = diagmul(dl, a);
    ComplexMatrix cl_expc = diagonalMatrix(dl, n1, n1) * a;
    
    std::cout << TAG_INFO << "Error (D_left * A - C_expc) = " << (cl - cl_expc).norm() << "\n";
    
    ComplexMatrix cr = diagmul(a, dr);
    ComplexMatrix cr_expc = a * diagonalMatrix(dr, n2, n2);
    
    std::cout << TAG_INFO << "Error (A * D_right - C_expc) = " << (cr - cr_expc).norm() << "\n";
    
    return 0;
}

/**
 * Test the solution of a linear system using LAPACK with multiple right-hand sides.
 */
int testLinearSolve() {
    std::cout << "====== TEST MATRIX SOLVE ======" << std::endl;
    
    const int n = 5;
    const int nrhs = 3;
    const double sigma = 1.;
    uint64_t seed = 1;
    
    ComplexMatrix a = gaussianRandomMatrix(n, n, sigma, seed);
    ComplexMatrix b = gaussianRandomMatrix(n, nrhs, sigma, seed+1);
    
    //a.print("Matrix A");
    //b.print("Independent term B");
    
    ComplexMatrix x(n, nrhs);
    
    solve(a, b, x);  // Call the linear solver (actually a wrapper for LAPACK's zgesv()).
    
    std::cout << TAG_INFO << "Error (A*X - B) = " << (a*x - b).norm() << "\n";
    
    return 0;
}

/**
 * Test the singular value decomposition using LAPACK.
 */
int testSVD() {
    std::cout << "====== TEST SVD ======" << std::endl;
    
    const int nrow = 10;
    const int ncol = 7;
    const int nsv = std::min(nrow, ncol);
    const double sigma = 1.;
    uint64_t seed = 1;
    
    ComplexMatrix a = gaussianRandomMatrix(nrow, ncol, sigma, seed);
    
    //a.print("Matrix A");
    
    ComplexMatrix u(nrow, nrow), vh(ncol, ncol);
    RealMatrix sdiag(nsv, 1);
    
    svd(a, sdiag, u, vh);  // Call the SVD function (actually a wrapper for LAPACK's zgesvd()).
    
    for (int i = 0; i < nsv-1; i++) {// First check that the singular values are ordered properly.
        if (sdiag(i,0) < sdiag(i+1,0)) {
            std::cout << TAG_WARN << "Singular value " << sdiag(i,0) << " is not sorted. Should be >=" << sdiag(i+1,0) << "." << std::endl;
        }
    }
    
    ComplexMatrix s = diagonalMatrix(sdiag, nrow, ncol);  // Construct a diagonal matrix from diagonal values.
    
    std::cout << TAG_INFO << "Error (U*S*V^H - A) = " << (u*s*vh - a).norm() << std::endl;
    std::cout << TAG_INFO << "Unitarity (U * U^H - 1) = " << (u*u.conj() - identityMatrix(nrow)).norm() << std::endl;
    std::cout << TAG_INFO << "Unitarity (V^H * V - 1) = " << (vh*vh.conj() - identityMatrix(ncol)).norm() << std::endl;
    
    return 0;
}

/**
 * Arbitrary function.
 */
dcomplex func(dcomplex z) {
    return std::sqrt(z);
}

/**
 * Test the matrix function application.
 */
int testApply() {
    std::cout << "====== TEST APPLY FUNCTION ======" << std::endl;
    
    const int n = 5;
    const double sigma = 1.;
    uint64_t seed = 1;
    
    ComplexMatrix a = gaussianRandomMatrix(n, n, sigma, seed);
    
    ComplexMatrix s = apply(&func, a);  // Apply a function to the matrix.
    
    std::cout << TAG_INFO << "Error (S*S - A) = " << (s*s - a).norm() << std::endl;
    
    return 0;
}

/**
 * Test the modal matrix.
 */
int testModalMatrix() {
    std::cout << "====== TEST MODAL MATRIX ======" << std::endl;
    
    const int n = 6;
    
    ComplexMatrix u = modalMatrix(n);
    
    u.print("Modal matrix U");
    
    std::cout << TAG_INFO << "Unitarity (U*U - 1) = " << (u*u - identityMatrix(n)).norm() << std::endl;
    
    return 0;
}

/**
 * Test the modal matrix.
 */
int testOpeningMatrix() {
    std::cout << "====== TEST OPENING MATRIX ======" << std::endl;
    
    const int n = 300;
    const double kh = PI/2.;
    
    auto start_build_old = std::chrono::steady_clock::now();
    
    ComplexMatrix a_old = openingMatrix_old(dcomplex(kh*kh, MEPS), n);
    
    double ctime_build_old = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_build_old).count();
    std::cout << TAG_INFO << "Build time (old) = " << ctime_build_old << " s.\n";
    
    auto start_build_new = std::chrono::steady_clock::now();
    
    ComplexMatrix a = openingMatrix(dcomplex(kh*kh, MEPS), n);
    
    double ctime_build_new = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_build_new).count();
    std::cout << TAG_INFO << "Build time (new) = " << ctime_build_new << " s. Speedup = " << ctime_build_old/ctime_build_new << ".\n";
    
    //a_old.print("Opening matrix A (old)");
    //a.print("Opening matrix A");
    
    std::cout << TAG_INFO << "Error (A - A_old) = " << (a - a_old).norm() << std::endl;
    
    return 0;
}

/**
 * Main function of the test.
 */
int main(int argc, char** argv) {
    
    //testHermitianConjugate();
    testMatrixProduct();
    testDiagmul();
    testLinearSolve();
    testSVD();
    testApply();
    
    testModalMatrix();
    testOpeningMatrix();
    
    return 0;
}
