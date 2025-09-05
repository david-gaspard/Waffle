/****
 * @date Created on 2025-08-05 at 15:58:49 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the implementation of the dense complex matrix.
 ***/
#include "ComplexMatrix.hpp"
#include "Image.hpp"
#include <random>
#include <iostream>

/**
 * Constructor of a matrix.
 */
ComplexMatrix::ComplexMatrix(const int nrow, const int ncol) {
    if (nrow <= 0) {
        throw std::invalid_argument("In ComplexMatrix(): Number of rows cannot be negative or zero.");
    }
    if (ncol <= 0) {
        throw std::invalid_argument("In ComplexMatrix(): Number of columns cannot be negative or zero.");
    }
    this->nrow = nrow;
    this->ncol = ncol;
    data = new dcomplex[nrow*ncol]();  // Allocate space for the matrix (column-major format) with zero initialization.
}

/**
 * Explicit copy constructor because of the underlying pointer "data".
 */
ComplexMatrix::ComplexMatrix(const ComplexMatrix& a) {
    nrow = a.nrow;
    ncol = a.ncol;
    const int nt = nrow*ncol;
    data = new dcomplex[nt];
    std::copy(a.data, a.data + nt, data);
}

/**
 * Conversion from real to complex matrix.
 */
ComplexMatrix::ComplexMatrix(const RealMatrix& a) {
    nrow = a.nrow;
    ncol = a.ncol;
    const int nt = nrow*ncol;
    data = new dcomplex[nt];
    std::copy(a.data, a.data + nt, data);  // Note that ComplexMatrix may access private data of RealMatrix.
}

/**
 * Destructor of a matrix. Frees the memory allocated by the constructor.
 */
ComplexMatrix::~ComplexMatrix() {
    delete[] data;
}

/***************************************************************************************************
 * GETTERS AND SETTERS
 ***************************************************************************************************/

/**
 * Returns the number of rows of this matrix.
 */
int ComplexMatrix::getNrow() const {
    return nrow;
}

/**
 * Returns the number of columns of this matrix.
 */
int ComplexMatrix::getNcol() const {
    return ncol;
}

/**
 * Returns a reference to the element (i, j) of the matrix.
 * The row index is "i", and the column index is "j".
 * Note that the indices "i" and "j" are based to zero (C/C++ convention).
 */
dcomplex& ComplexMatrix::operator()(const int i, const int j) const {
    if (i >= nrow || j >= ncol || i < 0 || j < 0) {
        std::string msg = "In ComplexMatrix(): Invalid index (i=" + std::to_string(i) + ", j=" + std::to_string(j) 
                        + ") given nrow=" + std::to_string(nrow) + " and ncol=" + std::to_string(ncol) + ".";
        throw std::out_of_range(msg);
    }
    return data[i + j*nrow];
}

/**
 * Overloads the assignement operator to perform a deep copy.
 */
ComplexMatrix& ComplexMatrix::operator=(const ComplexMatrix& a) {
    if (this == &a) {
        throw std::invalid_argument("In ComplexMatrix operator=(): Invalid self assignment.");
    }
    if (nrow != a.nrow || ncol != a.ncol) {
        std::string msg = "In ComplexMatrix operator=(): Matrices have different sizes. LHS = (" 
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
double ComplexMatrix::norm() const {
    
    double re, im, sumsq = 0.;
    
    for (int l = 0; l < nrow*ncol; l++) {// Loop over all elements.
        re = data[l].real();
        im = data[l].imag();
        sumsq += re*re + im*im;
    }
    
    return std::sqrt(sumsq);
}

/**
 * Returns the Hermitian conjugate (transpose conjugate) of the present matrix.
 */
ComplexMatrix ComplexMatrix::conj() const {
    
    ComplexMatrix c(ncol, nrow);
    
    for (int j = 0; j < ncol; j++) {// Loop over the columns.
        for (int i = 0; i < nrow; i++) {// Loop over the rows.
            c.data[j + i*ncol] = std::conj(data[i + j*nrow]);
        }
    }
    
    return c;
}

/**
 * Returns the real part of the present matrix.
 */
RealMatrix ComplexMatrix::real() const {
    
    RealMatrix re(nrow, ncol);
    
    for (int l = 0; l < nrow*ncol; l++) {// Loop over the matrix elements.
        re.data[l] = data[l].real();
    }
    
    return re;
}

/**
 * Returns the imaginary part of the present matrix.
 */
RealMatrix ComplexMatrix::imag() const {
    
    RealMatrix im(nrow, ncol);
    
    for (int l = 0; l < nrow*ncol; l++) {// Loop over the matrix elements.
        im.data[l] = data[l].imag();
    }
    
    return im;
}

/**
 * Defines the addition of two matrices of same dimensions.
 */
ComplexMatrix operator+(const ComplexMatrix& a, const ComplexMatrix& b) {
    
    if (a.nrow != b.nrow || a.ncol != b.ncol) {// Check for possible incompatibilities.
        std::string msg = "In operator+(): Incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    ComplexMatrix sum(a.nrow, a.ncol);
    
    for (int l = 0; l < a.nrow*a.ncol; l++) {
        sum.data[l] = a.data[l] + b.data[l];
    }
    
    return sum;
}

/**
 * Defines the subtraction of two matrices of same dimensions.
 */
ComplexMatrix operator-(const ComplexMatrix& a, const ComplexMatrix& b) {
    
    if (a.nrow != b.nrow || a.ncol != b.ncol) {// Check for possible incompatibilities.
        std::string msg = "In operator-(): Incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    ComplexMatrix diff(a.nrow, a.ncol);
    
    for (int l = 0; l < a.nrow*a.ncol; l++) {
        diff.data[l] = a.data[l] - b.data[l];
    }
    
    return diff;
}

/**
 * Multiplication of a matrix by a scalar.
 */
ComplexMatrix operator*(const dcomplex scalar, const ComplexMatrix& a) {
    
    ComplexMatrix b(a.nrow, a.ncol);
    
    for (int l = 0; l < a.nrow*a.ncol; l++) {// Loop on the matrix elements.
        b.data[l] = scalar * a.data[l];
    }
    
    return b;
}

/**
 * Multiplication of a matrix by a scalar the other way around.
 */
ComplexMatrix operator*(const ComplexMatrix& a, const dcomplex scalar) {
    return scalar*a;
}

/**
 * Product of complex matrices using BLAS.
 */
extern "C" void zgemm_(const char& transa, const char& transb, const int& m, const int& n, const int& k, const dcomplex& alpha, 
                       const dcomplex* a, const int& lda, const dcomplex* b, const int& ldb, const dcomplex& beta, dcomplex* c, const int& ldc);

/**
 * Defines the product of two matrices of compatible dimensions.
 * This function is a user-friendly wrapper to zgemm_().
 */
ComplexMatrix operator*(const ComplexMatrix& a, const ComplexMatrix& b) {
    
    if (a.ncol != b.nrow) {// Check for possible incompatibilities.
        std::string msg = "In operator*(): Product of incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    ComplexMatrix c(a.nrow, b.ncol);
    
    const char trans = 'N';   // Indicate that the matrices must not be transposed.
    const dcomplex alpha = dcomplex(1., 0.); // Constant "alpha" in C = alpha*A*B + beta*C (see "zgemm" in BLAS manual).
    const dcomplex beta = dcomplex(0., 0.);  // Constant "beta" in C = alpha*A*B + beta*C (see "zgemm" in BLAS manual).
    
    zgemm_(trans, trans, a.nrow, b.ncol, a.ncol, alpha, a.data, a.nrow, b.data, b.nrow, beta, c.data, a.nrow);
    
    return c;
}

/**
 * Multiplication by a diagonal matrix.
 * If A is a column matrix (A.ncol=1), then interpret A as a diagonal matrix and compute C_ij = A_i B_ij , 
 * otherwise if B is a row matrix (B.nrow=1), then interpret B as a diagonal matrix and compute C_ij = A_ij B_j .
 */
ComplexMatrix diagmul(const ComplexMatrix& a, const ComplexMatrix& b) {
    
    ComplexMatrix c(a.nrow, b.ncol);
    
    if (a.ncol == 1) {// Assume that A is the diagonal of a matrix multiplying B from the left.
        if (a.nrow != b.nrow) {
            std::string msg = "In diagmul(): Incompatible sizes of A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + "). Expected A.nrow=B.nrow.";
            throw std::invalid_argument(msg);
        }
        for (int l = 0; l < a.nrow*b.ncol; l++) {
            c.data[l] = a.data[l%a.nrow] * b.data[l];
        }
    }
    else if (b.nrow == 1) {// Assumes that B is the diagonal of a matrix multiplying A from the right.
        if (a.ncol != b.ncol) {
            std::string msg = "In diagmul(): Incompatible sizes of A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.nrow) + ", " + std::to_string(b.ncol) + "). Expected A.ncol=B.ncol.";
            throw std::invalid_argument(msg);
        }
        for (int l = 0; l < a.nrow*b.ncol; l++) {
            c.data[l] = a.data[l] * b.data[l/a.nrow];
        }
    }
    else {
        throw std::invalid_argument("In diagmul(): No diagonal in argument. Expected A.ncol=1 or B.nrow=1.");
    }
    return c;
}

/**
 * Solves a general complex linear system A.X = B using LAPACK.
 */
extern "C" void zgesv_(const int& n, const int& nrhs, dcomplex* a, const int& lda, int* ipiv, dcomplex* b, const int& ldb, int& info);
/**
 * Solves the linear system A.X = B for multiple right-hand sides using LAPACK.
 * This function is a user-friendly wrapper to zgesv_().
 */
void solve(const ComplexMatrix& a, const ComplexMatrix& b, ComplexMatrix& x) {
    
    // 1. First check for possible errors:
    if (a.nrow != a.ncol) {
        throw std::invalid_argument("In solve(): Cannot solve rectangular system.");
    }
    else if (a.ncol != x.nrow) {
        std::string msg = "In solve(): Invalid number of rows of X (unknown term), received x.nrow=" 
                        + std::to_string(x.nrow) + ", expected a.ncol=" + std::to_string(a.ncol) + ".";
        throw std::invalid_argument(msg);
    }
    else if (a.nrow != b.nrow) {
        std::string msg = "In solve(): Invalid number of rows of X (independent term), received b.nrow" 
                        + std::to_string(b.nrow) + ", expected a.nrow=" + std::to_string(a.nrow) + ".";
        throw std::invalid_argument(msg);
    }
    else if (b.ncol != x.ncol) {
        std::string msg = "In solve(): Incompatible number of columns of B and X, received b.ncol=" 
                        + std::to_string(b.ncol) + " and x.ncol=" + std::to_string(x.ncol) + ".";
        throw std::invalid_argument(msg);
    }
    
    // 2. Copy the matrix A (since it will be changed by LAPACK) and then call zgesv_():
    const int nta = a.nrow*a.ncol; // Total number of elements of A.
    const int ntb = b.nrow*b.ncol; // Total number of elements of B.
    
    auto ac = new dcomplex[nta];
    std::copy(a.data, a.data+nta, ac);
    std::copy(b.data, b.data+ntb, x.data);
    
    auto ipiv = new int[a.nrow];
    int info;
    
    zgesv_(a.nrow, b.ncol, ac, a.nrow, ipiv, x.data, b.nrow, info);  // Call the LAPACK solver.
    
    delete[] ac;  // Frees the allocated memory.
    delete[] ipiv;
    
    // 3. Check LAPACK's "info" for possible errors:
    if (info < 0) {
        std::string msg = "In solve(): Argument #" + std::to_string(-info) + " of zgesv_() had an illegal value.";
        throw std::runtime_error(msg);
    }
    else if (info > 0) {
        std::string msg = "In solve(): Element U(" + std::to_string(info) + ", " + std::to_string(info) + ") of LU factorization is exactly zero. The solution cound not be computed.";
        throw std::runtime_error(msg);
    }
}

/**
 * Compute the singular value decomposition A = U*S*V^H using LAPACK.
 */
extern "C" void zgesvd_(const char& jobu, const char& jobvt, const int& m, const int& n, dcomplex* a, const int& lda, double* s,
                        dcomplex* u, const int& ldu, dcomplex* vh, const int& ldvh, dcomplex* work, const int& lwork, double* rwork, int& info);
/**
 * Computes the singular value decomposition (SVD) of the given matrix A, i.e., the matrices U, S, and V such that A = U * S * V^H, 
 * where U and V are unitary matrices and S is a diagonal matrix with positive elements. Note that S is given as a vector (i.e., a Nx1 matrix).
 */
void svd(const ComplexMatrix& a, RealMatrix& s, ComplexMatrix& u, ComplexMatrix& vh) {
    
    // 1. Check for possible invalid arguments:
    const int nsv = std::min(a.nrow, a.ncol);  // Expected number of singular values.
    if (s.nrow != nsv || s.ncol != 1) {
        const std::string msg = "In svd(): Invalid dimensions of S, received (" + std::to_string(s.nrow) + ", " + std::to_string(s.ncol) 
                              + ") expected (" + std::to_string(nsv) + ", 1).";
        throw std::invalid_argument(msg);
    }
    else if (u.nrow != a.nrow || u.ncol != a.nrow) {
        const std::string msg = "In svd(): Invalid dimensions of U, received (" + std::to_string(u.nrow) + ", " + std::to_string(u.ncol) 
                              + ") expected (" + std::to_string(a.nrow) + ", " + std::to_string(a.nrow) + ").";
        throw std::invalid_argument(msg);
    }
    else if (vh.nrow != a.ncol || vh.ncol != a.ncol) {
        const std::string msg = "In svd(): Invalid dimensions of V^H, received (" + std::to_string(vh.nrow) + ", " + std::to_string(vh.ncol) 
                              + ") expected (" + std::to_string(a.ncol) + ", " + std::to_string(a.ncol) + ").";
        throw std::invalid_argument(msg);
    }
    
    // 2. Copy the matrix A (since it will be changed by LAPACK) and then call zgesvd_():
    const char job = 'A';  // Require the full SVD for both U and Vt.
    
    const int nta = a.nrow*a.ncol;  // Total number of matrix elements of A.
    auto ac = new dcomplex[nta];
    std::copy(a.data, a.data+nta, ac); // Copy the matrix A.
    //auto sreal = new double[nsv];
    const int lwork = std::max(2*nsv + std::max(a.nrow, a.ncol), 1);  // Minimum value of lwork (see LAPACK's zgesv() manual).
    auto work = new dcomplex[lwork];  // Allocate LAPACK's workspace.
    auto rwork = new double[5*nsv];
    int info;
    
    // Find the SVD using LAPACK:
    zgesvd_(job, job, a.nrow, a.ncol, ac, a.nrow, s.data, u.data, a.nrow, vh.data, a.ncol, work, lwork, rwork, info);
    
    //std::copy(sreal, sreal+nsv, s.data);
    
    delete[] ac;
    //delete[] sreal;
    delete[] work;
    delete[] rwork;
    
    // 3. Check LAPACK's "info" for possible errors:
    if (info < 0) {
        std::string msg = "In svd(): Argument #" + std::to_string(-info) + " of zgesvd_() had an illegal value.";
        throw std::runtime_error(msg);
    }
    else if (info > 0) {
        std::string msg = "In svd(): The QR factorization did not converge. The SVD could not be computed.";
        throw std::runtime_error(msg);
    }
}

/**
 * Compute the eigendecomposition A = U*W*U^-1 using LAPACK:
 */
extern "C" void zgeev_(const char& jobvl, const char& jobvr, const int& n, dcomplex* a, const int& lda, dcomplex* w,
                       dcomplex* vl, const int& ldvl, dcomplex* vr, const int& ldvr, dcomplex* work, const int& lwork, double* rwork, int& info);
/**
 * Apply a complex function f(z) on a complex matrix A using the eigendecomposition of this matrix, that is compute
 * f(A) = U*f(D)*U^-1 , where "U" is the matrix of right eigenvectors, and "D" is a diagonal matrix containing the eigenvalues.
 */
ComplexMatrix apply(dcomplex(*f)(dcomplex), const ComplexMatrix& a) {
    
    const int n = a.nrow;
    if (n != a.ncol) {// This function only makes sense for square matrices.
        throw std::invalid_argument("In apply(): Matrix must be square.");
    }
    
    // 1. First compute the eigenvalues/eigenvectors A*U = U*W :
    const char jobv = 'V', jobn = 'N';
    const int nt = n*n;  // Total number of matrix elements of A.
    ComplexMatrix fa(n, n);  // Final matrix will be used as temporary buffer.
    std::copy(a.data, a.data+nt, fa.data); // Copy the matrix A into A_tmp.
    auto u = new dcomplex[nt];  // Allocate space for the eigenvectors.
    auto w = new dcomplex[n];   // Allocate space for the eigenvalues.
    const int lwork = 2*n;
    auto work = new dcomplex[lwork];  // Allocate LAPACK's workspace (see LAPACK's zgeev() manual).
    auto rwork = new double[2*n];
    int info;
    
    zgeev_(jobn, jobv, n, fa.data, n, w, nullptr, n, u, n, work, lwork, rwork, info); // Compute the eigenvalues W and right eigenvectors U.
    
    if (info < 0) {
        std::string msg = "In apply(): Argument #" + std::to_string(-info) + " of zgeev_() had an illegal value.";
        throw std::runtime_error(msg);
    }
    else if (info > 0) {
        throw std::runtime_error("In apply(): The QR algorithm failed to compute all the eigenvalues. No eigenvectors have been computed.");
    }
    
    delete[] work;
    delete[] rwork;
    
    // 2. Then compute V = U^-1 :
    std::copy(u, u+nt, fa.data); // Copy the matrix U into A_tmp.
    auto ipiv = new int[n];
    auto v = new dcomplex[nt];
    for (int i = 0; i < n; i++) {// Set V as the identity matrix to compute V = U^-1.
        v[i*(n+1)] = dcomplex(1., 0.);
    }
    zgesv_(n, n, fa.data, n, ipiv, v, n, info); // Compute V = U^-1.
    
    delete[] ipiv;
    
    // 3. Compute V_new = f(W)*V :
    for (int l = 0; l < nt; l++) {
        v[l] = f(w[l%n]) * v[l];
    }
    
    // 4. Compute f(A) = U*V_new = U*f(W)*V :
    const dcomplex alpha = dcomplex(1., 0.); // Constant "alpha" in C = alpha*A*B + beta*C (see "zgemm" in BLAS manual).
    const dcomplex beta = dcomplex(0., 0.);  // Constant "beta" in C = alpha*A*B + beta*C (see "zgemm" in BLAS manual).
    
    zgemm_(jobn, jobn, n, n, n, alpha, u, n, v, n, beta, fa.data, n);
    
    delete[] u;
    delete[] v;
    delete[] w;
    
    return fa;
}

/***************************************************************************************************
 * CONSTANT MATRIX GENERATORS
 ***************************************************************************************************/

/**
 * Returns the identity matrix of given size.
 */
ComplexMatrix identityMatrix(const int n) {
    
    ComplexMatrix a(n, n);   // The constructor already checks for possible invalid inputs.
    
    for (int i = 0; i < n; i++) {// Loop over the nonzero elements.
        a(i, i) = dcomplex(1., 0.);
    }
    
    return a;
}

/**
 * Returns a diagonal matrix of dimensions (nrow, ncol) with the elements "diag" along the diagonal.
 */
ComplexMatrix diagonalMatrix(const ComplexMatrix& diag, const int nrow, const int ncol) {
    
    ComplexMatrix a(nrow, ncol);
    
    if (diag.getNcol() == 1) {// If A is a column vector, then interpret it as the diagonal.
        for (int i = 0; i < std::min(std::min(nrow, ncol), diag.getNrow()); i++) {
            a(i, i) = diag(i, 0);
        }
    }
    else if (diag.getNrow() == 1) {// If A is a row vector, then interpret it as the diagonal.
        for (int i = 0; i < std::min(std::min(nrow, ncol), diag.getNcol()); i++) {
            a(i, i) = diag(0, i);
        }
    }
    else {
        throw std::invalid_argument("In diagonalMatrix(): No diagonal in argument. Expected A.ncol=1 or A.nrow=1.");
    }
    return a;
}

/**
 * Returns a Gaussian random matrix of general complex kind with given dimensions and given standard deviation.
 * Argument "sigma" is the standard deviation of the real part (and imaginary part) of any of the matrix elements.
 */
ComplexMatrix gaussianRandomMatrix(const int nrow, const int ncol, const double sigma, const uint64_t seed) {
    
    if (sigma <= 0) {// Check for possible invalid arguments.
        throw std::invalid_argument("In gaussianRandomMatrix(): Standard deviation cannot be negative.");
    }
    
    ComplexMatrix a(nrow, ncol);   // The constructor already checks for possible invalid inputs.
    
    std::mt19937_64 rng;  // Instantiate the standard Mersenne Twister random number generator (64-bit-return version).
    rng.seed(seed);       // Initialize the random generator with the given seed.
    std::normal_distribution<double> random_normal(0., sigma);
    
    double x, y;
    
    for (int i = 0; i < nrow; i++) {// Loop over the rows.
        for (int j = 0; j < ncol; j++) {// Loop over the columns.
            x = random_normal(rng);  // This separation forces the random generator to be called in this specific order.
            y = random_normal(rng);
            a(i, j) = dcomplex(x, y);
        }
    }
    
    return a;
}

/**
 * Returns the discrete Laplacian matrix with -2 along the diagonal and +1 above and below the diagonal.
 * Note that the matrix assumes boundary conditions such that the function vanishes at the two ends.
 */
ComplexMatrix laplacianMatrix(const int n) {
    
    ComplexMatrix a(n, n);   // The constructor already checks for possible invalid inputs.
    
    a(0, 0) = dcomplex(-2., 0.);
    
    for (int i = 1; i < n; i++) {// Loop over the nonzero elements.
        a(i, i) = dcomplex(-2., 0.);
        a(i-1, i) = dcomplex(1., 0.);
        a(i, i-1) = dcomplex(1., 0.);
    }
    
    return a;
}

/**
 * Returns the i^th eigenvalue of the Laplacian matrix of size "n", which is just -4*sin((i + 1)*pi/(2*(n + 1)))^2.
 * Note that index "i" is zero based, i.e., it runs from 0 to n-1.
 */
double laplacianEigenvalue(const int i, const int n) {
    if (i < 0 || i >= n) {// Check for possible errors.
        std::string msg = "In laplacianEigenvalue(): Invalid index i=" + std::to_string(i) + ", expected in 0.." + std::to_string(n-1) + ".";
        throw std::invalid_argument(msg);
    }
    const double s = std::sin((i+1)*(PI/(2*(n+1))));
    return -4.*s*s;
}

/**
 * Returns the unitary matrix U containing the eigenvectors of the Laplacian matrix in columns.
 */
ComplexMatrix modalMatrix(const int n) {
    
    ComplexMatrix u(n, n);
    const double amp = std::sqrt(2./(n+1));
    double elem;
    
    for (int i = 0; i < n; i++) {// Loop over the rows.
        for (int j = 0; j <= i; j++) {// Loop over the columns in the lower triangle.
            elem = amp * std::sin((i+1)*(j+1)*(PI/(n+1)));
            u(i, j) = elem;   // Construct the matrix symmetrically.
            if (j != i) u(j, i) = elem;
        }
    }
    
    return u;
}

/**
 * Returns the matrix A involved in the openings (i.e., open boundary conditions), hence the name.
 * This matrix is exactly given by A = U * (-exp(-i*K_x)) * U^-1, where U is the modal matrix (containing the 
 * transverse waveguides modes in columns), and K_x is a diagonal matrix containing the effective wavenumbers of the lattice 
 * waveguide modes given by K_x,n = 2*arcsin(sqrt((k*h)^2 + d2ev_n)/2), where "d2ev" represents the eigenvalues of the 1D Laplacian matrix.
 * These eigenvalues are given by d2ev_n = -4*sin((n+1)*pi/(2*(N_y + 1)))^2. (see also function laplacianEigenvalues() above).
 * 
 * Arguments:
 * 
 * kh2  = Product kh*kh = (k*h)^2. Note that the imaginary part of "kh2" must be strictly positive (possibly vanishingly small).
 * n    = Size of the opening matrix.
 */
ComplexMatrix openingMatrix(const dcomplex kh2, const int n) {
    
    ComplexMatrix u = modalMatrix(n); // Generate the modal matrix. Time: O(N^2/2).
    ComplexMatrix opev(n, 1); // Computes the eigenvalues of the opening matrix.
    dcomplex d2ev;
    
    for (int i = 0; i < n; i++) {// Loop over the rows of "opev".
        d2ev = laplacianEigenvalue(i, n);  // Eigenvalue of the Laplacian matrix = -4*sin((i+1)*pi/(2*(n+1)))^2.
        // Eigenvalue of the opening matrix, -exp(-i*K_x) where K_x = 2*arcsin(sqrt(kh2 + d2ev)/2) :
        opev(i, 0) = -1. + (kh2 + d2ev)/2. + I*std::sqrt( (kh2 + d2ev) * (1. - (kh2 + d2ev)/4.) );
    }
    
    return u * diagmul(opev, u); // Multiplies the modal matrices together and returns. Time: O(N^3).
}

/**
 * Old version of the opening matrix generator.
 * 
 * @deprecated This version of the opening matrix generator is very slow. It is only used for testing.
 */
ComplexMatrix openingMatrix_old(const dcomplex kh2, const int n) {
    
    ComplexMatrix a(n, n);
    dcomplex elem, d2ev, opev;
    
    for (int i = 0; i < n; i++) {// Loop over the rows.
        for (int j = 0; j <= i; j++) {// Loop over the columns in the lower triangle.
            elem = 0.; // Compute the matrix element a(i, j).
            for (int k = 0; k < n; k++) {// Loop on the eigenvalues of the Laplacian matrix.
                d2ev = laplacianEigenvalue(k, n);  // Eigenvalue of the Laplacian matrix = -4*sin((k+1)*pi/(2*(n+1)))^2.
                // Eigenvalue of the opening matrix, -exp(-i*K_x) where K_x = 2*arcsin(sqrt(kh2 + d2ev)/2) :
                opev = -1. + (kh2 + d2ev)/2. + I*std::sqrt( (kh2 + d2ev) * (1. - (kh2 + d2ev)/4.) );
                elem += std::sin((i+1)*(k+1)*(PI/(n+1))) * std::sin((k+1)*(j+1)*(PI/(n+1))) * opev;
            }
            elem *= 2./(n+1);
            a(i, j) = elem;
            if (j != i) a(j, i) = elem;
        }
    }
    
    return a;
}

/***************************************************************************************************
 * PRINTING METHODS
 ***************************************************************************************************/

/**
 * Print the matrix to standard output.
 */
void ComplexMatrix::print(const std::string& name) const {
    
    dcomplex elem;
    
    std::cout << TAG_INFO << name << " = [\n";
    
    for (int i = 0; i < nrow; i++) {// Loop over the rows.
        std::cout << "\t[";
        for (int j = 0; j < ncol; j++) {// Loop over the columns.
            elem = data[i + j*nrow];
            std::cout << " " << elem.real() << (elem.imag() >= 0. ? "+" : "") << elem.imag() << "i ";
        }
        std::cout << "]\n";
    }
    std::cout << "]" << std::endl;
}

/**
 * Plots the complex matrix to a PNG file.
 */
void ComplexMatrix::savePNG(const std::string& filename) const {
    
    Image img(nrow, ncol);
    dcomplex elem;
    
    for (int i = 0; i < nrow; i++) {// Loop over the rows.
        for (int j = 0; j < ncol; j++) {// Loop over the columns.
            
            elem = data[i + j*nrow];  // Extract the matrix element (without changing the matrix state).
            
            //img(i, j) = complexColor1(elem);
            img(i, j) = complexColor2(elem);
        }
    }
    
    img.savePNG(filename);
}
