/****
 * @date Created on 2025-08-04 at 16:36:00 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the implementation for the SparseComplexMatrix object.
 ***/
#include "SparseComplexMatrix.hpp"
#include "Color.hpp"
#include "BaseTools.hpp"
#include <suitesparse/umfpack.h>
//#include <dmumps_c.h>
#include <algorithm>
#include <fstream>

/**
 * Constructor of the SparseComplexMatrix object when the dimensions are already known.
 */
SparseComplexMatrix::SparseComplexMatrix(const int nrow, const int ncol) {
    if (nrow <= 0) {
        throw std::invalid_argument("In setSize(): Number of rows cannot be negative.");
    }
    else if (ncol <= 0) {
        throw std::invalid_argument("In setSize(): Number of columns cannot be negative.");
    }
    this->nrow = nrow;
    this->ncol = ncol;
    symmetry = 0;
}

/**
 * Comparison function of two triplet. This function is used to sort the triplet vector of SparseMatrixObject in column-major ordering.
 * Returns "true" if triplet "t2" is located after "t1" according to the column-major ordering, "false" otherwise.
 */
bool compareTriplet(const Triplet& t1, const Triplet& t2) {
    return t1.j < t2.j || (t1.j == t2.j && t1.i < t2.i);
}

/**
 * Overload the streaming operator to print triplet automatically.
 */
std::ostream& operator<<(std::ostream& os, const Triplet& t) {
    os << "[ i=" << t.i << "  j=" << t.j << "  a=" << t.a.real() << (t.a.imag() >= 0. ? "+" : "") << t.a.imag() << "i ]";
    return os;
}

/*******************************************************************************
 * GETTERS AND SETTERS
 ******************************************************************************/

/**
 * Returns the number of nonzero elements in the SparseComplexMatrix.
 */
int64_t SparseComplexMatrix::getNnz() const {
    return triplet.size();
}

/**
 * Returns the number of rows of the sparse matrix.
 */
int SparseComplexMatrix::getNrow() const {
    return nrow;
}

/**
 * Returns the number of columns of the sparse matrix.
 */
int SparseComplexMatrix::getNcol() const {
    return ncol;
}

/**
 * Returns the density of the present sparse matrix normalized to 1.
 * The maximum density is 1.
 */
double SparseComplexMatrix::density() const {
    return ((double)triplet.size())/(nrow*ncol);
}

 /**
  * Returns the symmetry flag of the sparse matrix. 0=Nonsymmetric, 1=Symmetric.
  */   
int SparseComplexMatrix::getSymmetry() const {
    return symmetry;
}

/**
 * ASsigns the flag corresponding to the symmetry index.
 */
void SparseComplexMatrix::setSymmetry(const int symmetry) {
    this->symmetry = symmetry;
}

/**
 * Check that the given indices are valid.
 */
void SparseComplexMatrix::checkIndices(const int i, const int j) const {
    
    if (i >= nrow || i < 0) {
        std::string msg = "In SparseComplexMatrix(): Invalid index 'i' in (i=" + std::to_string(i) + ", j=" + std::to_string(j) 
                        + ") given the matrix size (" + std::to_string(nrow) + ", " + std::to_string(ncol) 
                        + "). Expected 'i' in range 0.." + std::to_string(nrow-1) + ".";
        throw std::out_of_range(msg);
    }
    else if (j >= ncol || j < 0) {
        std::string msg = "In SparseComplexMatrix(): Invalid index 'j' in (i=" + std::to_string(i) + ", j=" + std::to_string(j) 
                        + ") given the matrix size (" + std::to_string(nrow) + ", " + std::to_string(ncol) 
                        + "). Expected 'j' in range 0.." + std::to_string(ncol-1) + ".";
        throw std::out_of_range(msg);
    }
}

/**
 * Returns the element (i, j) of the matrix, or zero if the element does not exist.
 * This method is a getter: It does not change the matrix state.
 */
dcomplex SparseComplexMatrix::get(const int i, const int j) const {
    
    checkIndices(i, j);
    
    Triplet t;
    t.i = i; t.j = j;
    
    auto p = std::lower_bound(triplet.begin(), triplet.end(), t, compareTriplet);  // Find the first element >= than the given triplet using binary search.
    
    if (p == triplet.end() || p->i != i || p->j != j) {// If the triplet does not exist in the current matrix.
        // Note that the first condition prevents the pointer "p" to be dereferenced out of range (because end() always points out of the vector).
        return dcomplex(0., 0.);
    }
    
    return p->a;
}

/**
 * Returns the address of the matrix element A_ij. If the matrix element does not exist, then creates
 * a new matrix element before returning the address of the matrix element.
 * The triplets are automatically sorted so that binary search can be used to reduce the access time to O(log2(nnz)).
 */
dcomplex& SparseComplexMatrix::operator()(const int i, const int j) {
    
    checkIndices(i, j);
    
    Triplet t;
    t.i = i; t.j = j; t.a = dcomplex(0., 0.);
    
    auto p = std::lower_bound(triplet.begin(), triplet.end(), t, compareTriplet);  // Find the first element >= than the given triplet using binary search.
    
    if (p == triplet.end() || p->i != i || p->j != j) {// If the triplet does not exist in the current matrix.
        // Note that the first condition prevents the pointer "p" to be dereferenced out of range (because end() always points out of the vector).
        //std::cout << TAG_INFO << "Inserting new triplet " << t << "." << std::endl;
        p = triplet.insert(p, t);  // Insert the new element. Pointer "p" is now guaranteed to point to the added triplet.
    }
    
    return p->a;
}

/*******************************************************************************
 * PRINTING METHODS
 ******************************************************************************/

/**
 * Print essential information about the sparse matrix without showing the full content.
 */
void SparseComplexMatrix::printInfo(const std::string& name) const {
    std::cout << TAG_INFO << "Sparse complex matrix '" << name << "': Nrow=" << nrow << ", Ncol=" << ncol << ", Nnz=" << triplet.size() 
              << ", Density=" << ((100.*triplet.size())/(nrow*ncol)) << "%, Symmetric=" << (isSymmetric() ? "true" : "false") << ".\n";
}

/**
 * Prints the vector of triplets to std output.
 */
void SparseComplexMatrix::print(const std::string& name) const {
    std::cout << TAG_INFO << name << " (" << nrow << "x" << ncol << ", nnz=" << triplet.size() 
              << ", density=" << ((100.*triplet.size())/(nrow*ncol)) << "%) = [\n";
    for (Triplet t : triplet) {// Loop on the triplets.
        std::cout << "\t" << t << "\n";
    }
    std::cout << "]" << std::endl;
}

/**
 * Prints the sparsity pattern of the present sparse matrix to a portable pixmap file, a PPM file (see: https://en.wikipedia.org/wiki/Netpbm).
 */
void SparseComplexMatrix::printImage(const std::string& filename) const {
    
    // Open a file in binary mode:
    std::ofstream ofs(filename, std::ios::binary);
    ofs << "P6 " << (uint32_t)ncol << ' ' << (uint32_t)nrow << ' ' << (uint32_t)MAX_COLOR << ' ';
    
    const std::string msg = "Image sparse matrix";
    auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    dcomplex elem;
    for (int i = 0; i < nrow; i++) {// Loop over the rows.
        for (int j = 0; j < ncol; j++) {// Loop over the columns.
            
            elem = get(i, j);  // Extract the matrix element (without changing the matrix state).
            
            if (elem.real() != 0. || elem.imag() != 0.) {// If the matrix element is non-zero, then pick up a color.
                ofs << complexColor1(elem);
                //ofs << complexColor2(elem);
            }
            else {// Otherwise, draw white pixel.
                ofs << COLOR_WHITE;
            }
        }
        if (i % (nrow/200) == 0) {// Display a progress bar because this may take a while.
            printProgressBar(i+1, nrow, msg, start);
        }
    }
    endProgressBar(start); // Returns the computation time.
}

/*******************************************************************************
 * MATHEMATICAL OPERATIONS
 ******************************************************************************/

/**
 * Computes the Frobenius norm of the current sparse matrix.
 */
double SparseComplexMatrix::norm() const {
    
    double re, im, sumsq = 0.;
    
    for (Triplet t : triplet) {// Loop over the triplets.
        re = t.a.real();
        im = t.a.imag();
        sumsq += re*re + im*im;
    }
    
    return std::sqrt(sumsq);
}

/**
 * Returns the Hermitian conjugate (i.e., conjugate transpose) of the current sparse matrix.
 */
SparseComplexMatrix SparseComplexMatrix::conj() const {
    
    SparseComplexMatrix c(ncol, nrow);
    
    for (Triplet t : triplet) {// Loop over the triplets.
        c(t.j, t.i) = std::conj(t.a);
    }
    
    return c;
}

/**
 * Returns "true" is the matrix is symmetric, "false" otherwise.
 */
bool SparseComplexMatrix::isSymmetric() const {
    
    if (nrow != ncol) {// If the matrix is non-square, then skip further calculations.
        return false;
    }
    
    Triplet test;
    
    for (Triplet t : triplet) {// Loop over the triplets.
        if (t.i > t.j) {// If the element is in the lower triangle.
            test.i = t.j; test.j = t.i; test.a = t.a;  // Initialize triplet "test" with reversed indices.
            auto p = std::lower_bound(triplet.begin(), triplet.end(), test, compareTriplet);  // Find the first element >= than "test" using binary search.
            
            if (p == triplet.end() || p->i != test.i || p->j != test.j || p->a != test.a) {
                return false;
            }
        }
    }
    
    return true;
}

/**
 * Multiplication of a sparse matrix by a scalar factor.
 */
SparseComplexMatrix operator*(const dcomplex scalar, const SparseComplexMatrix& a) {
    
    SparseComplexMatrix b(a.nrow, a.ncol);
    
    for (Triplet t : a.triplet) {// Loop on the triplets of A.
        b(t.i, t.j) = scalar * t.a;
    }
    
    return b;
}

/**
 * Multiplication of a sparse matrix by a scalar the other way around.
 */
SparseComplexMatrix operator*(const SparseComplexMatrix& a, const dcomplex scalar) {
    return scalar*a;
}

/**
 * Multiplication of a sparse matrix by a dense matrix
 */
ComplexMatrix operator*(const SparseComplexMatrix& a, const ComplexMatrix& b) {
    
    if (a.ncol != b.getNrow()) {// Check for possible incompatibilities.
        std::string msg = "In operator*(): Product of incompatible matrix sizes A(" 
                        + std::to_string(a.nrow) + ", " + std::to_string(a.ncol) + ") and B(" 
                        + std::to_string(b.getNrow()) + ", " + std::to_string(b.getNcol()) + ").";
        throw std::invalid_argument(msg);
    }
    
    ComplexMatrix c(a.nrow, b.getNcol());  // Declare a dense matrix (initially zero).
    
    for (Triplet t : a.triplet) {// Loop over the triplets of A.
        
        for (int k = 0; k < b.getNcol(); k++) {// Loop over the columns of B.
            c(t.i, k) += t.a * b(t.j, k);   // C_ik = Sum_j A_ij B_jk
        }
        
    }
    
    return c;
}

/**
 * Check if the dimensions of the given matrices are compatible in view of solving the linear system A.x = b.
 */
void checkSolverInput(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x) {
    if (a.getNrow() != a.getNcol()) {
        std::string msg = "In checkSolverInput(): Matrix A must be square, received A.nrow=" 
                        + std::to_string(a.getNrow()) + " != A.ncol=" + std::to_string(a.getNcol()) + ".";
        throw std::invalid_argument(msg);
    }
    else if (a.getNcol() != x.getNrow()) {
        std::string msg = "In checkSolverInput(): Invalid number of rows of X (unknown term), received X.nrow=" 
                        + std::to_string(x.getNrow()) + ", expected A.ncol=" + std::to_string(a.getNcol()) + ".";
        throw std::invalid_argument(msg);
    }
    else if (a.getNrow() != b.getNrow()) {
        std::string msg = "In checkSolverInput(): Invalid number of rows of B (independent term), received B.nrow" 
                        + std::to_string(b.getNrow()) + ", expected A.nrow=" + std::to_string(a.getNrow()) + ".";
        throw std::invalid_argument(msg);
    }
    else if (b.getNcol() != x.getNcol()) {
        std::string msg = "In checkSolverInput(): Incompatible number of columns of B and X, received B.ncol=" 
                        + std::to_string(b.getNcol()) + " and X.ncol=" + std::to_string(x.getNcol()) + ".";
        throw std::invalid_argument(msg);
    }
}

/**
 * Solve the sparse linear system A.X = B for multiple right-hand sides using the UMFPACK solver.
 * 
 * Arguments:
 * 
 * A  = Sparse complex matrix object (unchanged).
 * B  = Sparse complex matrix containing the right-hand sides (unchanged).
 * X  = Complex matrix containing the solutions on output.
 */
void solveUmfpack(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x) {
    
    // 1. First check for possible errors in the input:
    checkSolverInput(a, b, x);
    
    // 2. Convert the matrix into UMFPACK's compressed column format:
    int32_t nnz = a.triplet.size();  // Estimate the number of nonzero (with/without symmetry).
    auto pcol = new int32_t[a.ncol+1]();
    auto irow = new int32_t[nnz];
    auto real = new double[nnz];
    auto imag = new double[nnz]; 
    
    nnz = 0;  // Reset the number of nonzero to count them.
    int32_t jcol = 0;  // Current column index.
    
    for (Triplet t : a.triplet) {// Loop on the triplets in column-major order.
        
        if (t.j != jcol) {
            jcol = t.j;
            pcol[jcol] = nnz;
        }
        irow[nnz] = t.i;
        real[nnz] = t.a.real();
        imag[nnz] = t.a.imag();
        nnz++;
    }
    pcol[a.ncol] = nnz;
    
    // 4. Call UMFPACK for each right-hand side:
    double info[UMFPACK_INFO];  // Output information (including UMFPACK's returned value "status").
    double control[UMFPACK_CONTROL]; // UMFPACK's control parameters (input).
    umfpack_zi_defaults(control); // Setup the default UMFPACK parameters for complex arrays (and long indices).
    control[UMFPACK_PRL] = 1;     // Control UMFPACK's default printing level (see UMFPACK's Manual).
    control[UMFPACK_IRSTEP] = 0;  // Set number of iterative refinement to zero. 
    
    //umfpack_zi_report_matrix(a.n, a.n, pcol, irow, real, imag, 1, control);
    
    auto worki = new int32_t[a.ncol];  // Workspace used by UMFPACK (move to class field).
    auto work = new double[4*a.ncol];  // Workspace array. Size=4*ncol without iterative refinement, and size=10*ncol with iterative refinement (see manual).
    
    void *symbolic, *numeric;  // Symbolic and numeric factorization of the sparse matrix.
    
    umfpack_zi_symbolic(a.nrow, a.ncol, pcol, irow, real, imag, &symbolic, control, info);  // Perform symbolic reording to minimize fill-in.
    umfpack_zi_numeric(pcol, irow, real, imag, symbolic, &numeric, control, info); // Perform numericl LU factorization.
    umfpack_zi_free_symbolic(&symbolic); // Free the memory allocated by the symbolic factorization.
    
    auto breal = new double[b.getNrow()];
    auto bimag = new double[b.getNrow()];
    auto xreal = new double[x.getNrow()];
    auto ximag = new double[x.getNrow()];
    
    for (int jrhs = 0; jrhs < b.getNcol(); jrhs++) {// Loop over the right-hand sides.
        
        for (int i = 0; i < b.getNrow(); i++) {// Copy the right-hand side in "breal" and "bimag".
            breal[i] = b.get(i, jrhs).real();
            bimag[i] = b.get(i, jrhs).imag();
        }
        
        umfpack_zi_wsolve(UMFPACK_A, pcol, irow, real, imag, xreal, ximag, breal, bimag, numeric, control, info, worki, work);
        
        for (int i = 0; i < x.getNrow(); i++) {// Copy the solution in "x".
            x(i, jrhs) = dcomplex(xreal[i], ximag[i]);
        }
    }
    
    umfpack_zi_free_numeric(&numeric); // Free the memory allocated by the numeric factorization.
    
    // Free the memory:
    delete[] irow;
    delete[] pcol;
    delete[] real;
    delete[] imag;
    delete[] worki;
    delete[] work;
    delete[] xreal;
    delete[] ximag;
    delete[] breal;
    delete[] bimag;
}

/**
 * Solve the sparse linear system A.X = B for multiple right-hand sides using the UMFPACK solver.
 * 
 * @deprecated Old version using too much memory.............
 * 
 * Arguments:
 * 
 * A  = Sparse complex matrix object (unchanged).
 * B  = Complex matrix containing the right-hand sides (unchanged).
 * X  = Complex matrix containing the solutions on output.
 */
void solveUmfpack_old(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x) {
    
    // 1. First check for possible errors in the input:
    checkSolverInput(a, b, x);
    
    // 2. Convert the matrix into UMFPACK's triplet format:
    int32_t nnz = a.triplet.size();  // Estimate the number of nonzero (with/without symmetry).
    if (a.symmetry != 0) {// If the matrix is declared as symmetric, then double the spacing.
        nnz *= 2;
    }
    auto ti = new int32_t[nnz];
    auto tj = new int32_t[nnz];
    auto treal = new double[nnz];
    auto timag = new double[nnz];
    
    nnz = 0;  // Reset the number of nonzero to count them.
    
    for (Triplet t : a.triplet) {// Loop on the triplets.
        
        ti[nnz] = t.i;  // Subtract the base index because UMFPACK uses zero-based indices.
        tj[nnz] = t.j;
        treal[nnz] = t.a.real();
        timag[nnz] = t.a.imag();
        nnz++;
        
        if (a.symmetry != 0 && t.i != t.j) {// If the matrix is symmetric, then add the symmetric elements.
            ti[nnz] = t.j;  // Note the reversal of indices.
            tj[nnz] = t.i;
            treal[nnz] = t.a.real();
            timag[nnz] = t.a.imag();
            nnz++;
        }
    }
    
    // 3. Convert the matrix into UMFPACK's compressed column format:
    auto irow = new int32_t[nnz];
    auto pcol = new int32_t[a.ncol+1];
    auto real = new double[nnz];
    auto imag = new double[nnz]; 
    
    int status = umfpack_zi_triplet_to_col(a.nrow, a.ncol, nnz, ti, tj, treal, timag, pcol, irow, real, imag, nullptr);
    
    if (status != UMFPACK_OK) {// Check for possible UMFPACK errors.
        std::string info = "In solveUmfpack(): UMFPACK's triplet_to_col() failed. Error code = " + std::to_string(status) + ".";
        throw std::runtime_error(info);
    }
    
    delete[] ti;
    delete[] tj;
    delete[] treal;
    delete[] timag;
    
    // 4. Call UMFPACK for each right-hand side:
    double info[UMFPACK_INFO];  // Output information (including UMFPACK's returned value "status").
    double control[UMFPACK_CONTROL]; // UMFPACK's control parameters (input).
    umfpack_zi_defaults(control); // Setup the default UMFPACK parameters for complex arrays (and long indices).
    control[UMFPACK_PRL] = 5;     // Control UMFPACK's default printing level (see UMFPACK's Manual).
    control[UMFPACK_IRSTEP] = 0;  // Set number of iterative refinement to zero. 
    
    //umfpack_zi_report_matrix(a.n, a.n, pcol, irow, real, imag, 1, control);
    
    auto worki = new int32_t[a.ncol];  // Workspace used by UMFPACK (move to class field).
    auto work = new double[4*a.ncol];  // Workspace array. Size=4*ncol without iterative refinement, and size=10*ncol with iterative refinement (see manual).
    
    void *symbolic, *numeric;  // Symbolic and numeric factorization of the sparse matrix.
    
    umfpack_zi_symbolic(a.nrow, a.ncol, pcol, irow, real, imag, &symbolic, control, info);  // Perform symbolic reording to minimize fill-in.
    umfpack_zi_numeric(pcol, irow, real, imag, symbolic, &numeric, control, info); // Perform numericl LU factorization.
    umfpack_zi_free_symbolic(&symbolic); // Free the memory allocated by the symbolic factorization.
    
    auto breal = new double[b.getNrow()];
    auto bimag = new double[b.getNrow()];
    auto xreal = new double[x.getNrow()];
    auto ximag = new double[x.getNrow()];
    
    for (int jrhs = 0; jrhs < b.getNcol(); jrhs++) {// Loop over the right-hand sides.
        
        for (int i = 0; i < b.getNrow(); i++) {// Copy the right-hand side in "breal" and "bimag".
            breal[i] = b.get(i, jrhs).real();
            bimag[i] = b.get(i, jrhs).imag();
        }
        
        umfpack_zi_wsolve(UMFPACK_A, pcol, irow, real, imag, xreal, ximag, breal, bimag, numeric, control, info, worki, work);
        
        for (int i = 0; i < x.getNrow(); i++) {// Copy the solution in "x".
            x(i, jrhs) = dcomplex(xreal[i], ximag[i]);
        }
    }
    
    umfpack_zi_free_numeric(&numeric); // Free the memory allocated by the numeric factorization.
    
    // Free the memory:
    delete[] irow;
    delete[] pcol;
    delete[] real;
    delete[] imag;
    delete[] worki;
    delete[] work;
    delete[] xreal;
    delete[] ximag;
    delete[] breal;
    delete[] bimag;
}

/**
 * Solve the sparse linear system A.X = B for multiple right-hand sides using the MUMPS solver.
 * 
 * Arguments:
 * 
 * A  = Sparse complex matrix object (unchanged).
 * B  = Sparse complex matrix containing the right-hand sides (unchanged).
 * X  = Complex matrix containing the solutions on output.
 */
void solveMumps(const SparseComplexMatrix& a, const SparseComplexMatrix& b, ComplexMatrix& x) {
    
    // 1. First check for possible errors in the input:
    checkSolverInput(a, b, x);
    
    std::cout << TAG_WARN << "In solveMumps(): Not implemented yet !\n";
    
    // TODO: See p.121 of MUMPS manual: https://mumps-solver.org/doc/userguide_5.8.1.pdf  ......................
    
}
