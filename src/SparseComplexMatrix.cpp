/****
 * @date Created on 2025-08-27 at 14:58:13 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the implementation for the SparseComplexMatrix object.
 ***/
#include "SparseComplexMatrix.hpp"
#include "Color.hpp"
#include "BaseTools.hpp"
#include <suitesparse/umfpack.h>
//#include <dmumps_c.h>
#include <iostream>

/**
 * Constructor of the SparseComplexMatrix object when the dimensions are already known.
 */
SparseComplexMatrix::SparseComplexMatrix(const int nrow, const int ncol) {
    if (nrow <= 0) {
        throw std::invalid_argument("In SparseComplexMatrix(): Number of rows cannot be negative.");
    }
    else if (ncol <= 0) {
        throw std::invalid_argument("In SparseComplexMatrix(): Number of columns cannot be negative.");
    }
    this->nrow = nrow;
    this->ncol = ncol;
}

/**
 * Overloads the comparison operator between two indices. This overload is necessary to sort the matrix elements of the Map in column-major ordering.
 * Returns "true" if indices "id2" is located after "id1" according to the column-major ordering, "false" otherwise.
 */
bool operator<(const Indices& id1, const Indices& id2) {
    return id1.j < id2.j || (id1.j == id2.j && id1.i < id2.i);
}

/*******************************************************************************
 * GETTERS AND SETTERS
 ******************************************************************************/

/**
 * Returns the number of nonzero elements in the SparseComplexMatrix.
 */
int64_t SparseComplexMatrix::getNnz() const {
    return data.size();
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
 * Returns the density of the present sparse matrix normalized to 1. The maximum density is 1.
 * The static cast is needed because the product nrow*ncol can exceed the capability of int32.
 */
double SparseComplexMatrix::density() const {
    return static_cast<double>(data.size())/(static_cast<double>(nrow) * static_cast<double>(ncol));
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
 * Returns the address of the matrix element A_ij. If the matrix element does not exist, then creates
 * a new matrix element before returning the address of the matrix element.
 * The triplets are automatically sorted so that binary search can be used to reduce the access time to O(log2(nnz)).
 */
dcomplex& SparseComplexMatrix::operator()(const int i, const int j) {
    
    checkIndices(i, j);
    
    return data[Indices{i, j}];
}

/**
 * Returns the element (i, j) of the matrix, or zero if the element does not exist.
 * This method is a getter: In contrast to operator(), it does not change the matrix state.
 */
dcomplex SparseComplexMatrix::get(const int i, const int j) const {
    
    checkIndices(i, j);
    
    auto p = data.find(Indices{i, j});
    
    if (p == data.end()) {// If the matrix element is not found, then return zero.
        return 0.;
    }
    
    return p->second;  // If the element exists, then return it.
}

/*******************************************************************************
 * PRINTING METHODS
 ******************************************************************************/

/**
 * Print essential information about the sparse matrix without showing the full content.
 */
void SparseComplexMatrix::summary(const std::string& name) const {
    std::cout << TAG_INFO << "Sparse complex matrix '" << name << "': Nrow=" << nrow << ", Ncol=" << ncol << ", Nnz=" << data.size() 
              << ", Density=" << density() << "%, Symmetric=" << (isSymmetric() ? "true" : "false") << ".\n";
}

/**
 * Prints all the matrix elements to std output.
 */
void SparseComplexMatrix::print(const std::string& name) const {
    std::cout << TAG_INFO << name << " (" << nrow << "x" << ncol << ", nnz=" << data.size() 
              << ", density=" << density() << "%) = [\n";
    for (const auto& [id, elem] : data) {// Loop on the nonzero matrix elements.
        std::cout << "\t(" << id.i << ", " << id.j << ") = " << elem.real() << (elem.imag() >= 0. ? "+" : "") << elem.imag() << "i\n";
    }
    std::cout << "]" << std::endl;
}

/**
 * Prints the sparsity pattern of the present sparse matrix to a portable pixmap file, a PPM file (see: https://en.wikipedia.org/wiki/Netpbm).
 */
void SparseComplexMatrix::saveImage(const std::string& filename) const {
    
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
    
    for (const auto& [id, elem] : data) {// Loop over the nonzero matrix elements.
        re = elem.real();
        im = elem.imag();
        sumsq += re*re + im*im;
    }
    
    return std::sqrt(sumsq);
}

/**
 * Returns the Hermitian conjugate (i.e., conjugate transpose) of the current sparse matrix.
 */
SparseComplexMatrix SparseComplexMatrix::conj() const {
    
    SparseComplexMatrix c(ncol, nrow);
    
    for (const auto& [id, elem] : data) {// Loop over the matrix elements.
        c(id.j, id.i) = std::conj(elem);
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
    
    for (const auto& [id, elem] : data) {// Loop over the nonzero matrix elements.
        if (id.i > id.j) {// If the element is in the lower triangle.
            
            auto p = data.find(Indices{id.j, id.i});  // Find the matrix element with permuted indices (j, i).
            
            if (p == data.end() || p->second != elem) {// If the symmetric element does not exist or is different, then return "false".
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
    
    SparseComplexMatrix b(a); // Copy constructor.
    
    for (auto& [id, elem] : b.data) {// Loop over the nonzero matrix elements of A.
        //b(id.i, id.j) = scalar * elem;
        elem *= scalar;
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
    
    ComplexMatrix c(a.nrow, b.getNcol());  // Declare a dense matrix C = A*B (initially zero).
    
    for (const auto& [id, elem] : a.data) {// Loop over the triplets of A.
        
        for (int k = 0; k < b.getNcol(); k++) {// Loop over the columns of B.
            c(id.i, k) += elem * b(id.j, k);   // C_ik = Sum_j A_ij B_jk
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
    int32_t nnz = a.data.size();  // Estimate the number of nonzero (with/without symmetry).
    auto pcol = new int32_t[a.ncol+1]();
    auto irow = new int32_t[nnz];
    auto real = new double[nnz];
    auto imag = new double[nnz]; 
    
    nnz = 0;  // Reset the number of nonzero to count them.
    int32_t jcol = 0;  // Current column index.
    
    for (const auto& [id, elem] : a.data) {// Loop on the triplets in column-major order.
        
        if (id.j != jcol) {
            jcol = id.j;
            pcol[jcol] = nnz;
        }
        irow[nnz] = id.i;
        real[nnz] = elem.real();
        imag[nnz] = elem.imag();
        nnz++;
    }
    pcol[a.ncol] = nnz;
    
    // 4. Call UMFPACK for each right-hand side:
    double info[UMFPACK_INFO];  // Output information (including UMFPACK's returned value "status").
    double control[UMFPACK_CONTROL]; // UMFPACK's control parameters (input).
    umfpack_zi_defaults(control); // Setup the default UMFPACK parameters for complex arrays (and long indices).
    control[UMFPACK_PRL] = 1;     // Control UMFPACK's default printing level (see UMFPACK's Manual).
    control[UMFPACK_IRSTEP] = 0;  // Set number of iterative refinement to zero. Gives sensible acceleration.
    
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
