/****
 * @date Created on 2025-08-29 at 13:23:26 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the implementation for the SparseComplexMatrix object.
 ***/
#include "SparseComplexMatrix.hpp"
#include "Image.hpp"
#include "BaseTools.hpp"
#include <suitesparse/umfpack.h>
#include <zmumps_c.h>  // MUMPS header for double precision complex arithmetic (code tested with MUMPS 5.8.0).
#include <algorithm>
#include <fstream>

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
    sorted = false;  //At the beginning, the triplets are not sorted.
}

/**
 * Overloads the comparison operator between two triplets. This overload is necessary to sort the matrix elements in column-major ordering.
 * Returns "true" if triplet "t2" is located after "t1" according to the column-major ordering, "false" otherwise.
 */
bool operator<(const Triplet& t1, const Triplet& t2) {
    return t1.j < t2.j || (t1.j == t2.j && t1.i < t2.i);
}

/**
 * Overloads the equality operator between two triplets. This overload is necessary to remove the duplicate matrix elements (see: finalize()).
 * Returns "true" if triplet "t1" is located at the same position as "t1", "false" otherwise.
 */
bool operator==(const Triplet& t1, const Triplet& t2) {
    return t1.i == t2.i && t1.j == t2.j;
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
    return static_cast<double>(triplet.size())/(static_cast<double>(nrow) * static_cast<double>(ncol));
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
 * Check if the present matrix if sorted, throws an error if it is not.
 */
void SparseComplexMatrix::checkSorted(const std::string& name) const {
    if (not sorted) {
        throw std::logic_error("In " + name + ": SparseComplexMatrix is not completely initialized. Please use finalize().");
    }
}

/**
 * Preallocate the memory occupied by the nonzero matrix elements.
 * This operation is not necessary but recommended before adding matrix elements because it avoids the cost of reallocation.
 */
void SparseComplexMatrix::allocate(const int nnz) {
    if (nnz <= 0) {
        std::string msg = "In allocate(): Sparse matrix allocation size must be positive (received " + std::to_string(nnz) + ").";
        throw std::invalid_argument(msg);
    }
    triplet.reserve(nnz);
}

/**
 * Returns the address of the matrix element A_ij. If the matrix element does not exist, then creates
 * a new matrix element before returning the address of the matrix element.
 * The method is implemented in such a way that, if the elements are NORT sorted (in building phase), elements are pushed back (O(1) in time),
 * otherwise if the elements are sorted (after finalize), binary search is used (which is O(log(N)) in time).
 * In the exceptional case where an element is added in an already sorted matrix, insertion is used to maintain ordering (but it is O(N) in time).
 */
dcomplex& SparseComplexMatrix::operator()(const int i, const int j) {
    
    checkIndices(i, j);  // First check the validity of indices.
    // DO NOT call finalize() here because it will degrade the construction speed considerably, passing from O(N) to O(N^2).
    
    Triplet t;
    t.i = i; t.j = j; t.a = dcomplex(0., 0.);
    
    if (not sorted) {// If the matrix is not sorted, i.e., still in construction, then use push_back(), it is O(1) in time.
        triplet.push_back(t);  // Using push_back() is O(1) in time (and very fast), but O(N) if reallocation happens.
                               // To mitigate this, just prealloce the vector "triplet" with the known number of nonzero elements.
        return triplet.back().a;  // Returns the reference to the added matrix element.
    }
    else {// If the matrix is sorted, then uses binary search which is O(N*log(N)) in time.
        auto p = std::lower_bound(triplet.begin(), triplet.end(), t); // Find the first element >= than the given triplet using binary search.
        
        if (p == triplet.end() || p->i != i || p->j != j) {// If the triplet does not exist in the current matrix.
            // Note that the first condition prevents the pointer "p" to be dereferenced out of range (because end() always points out of the vector).
            //std::cout << TAG_INFO << "Inserting new triplet " << t << "." << std::endl;
            p = triplet.insert(p, t);  // Insert the new element. Pointer "p" is now guaranteed to point to the added triplet. This is O(N) in time.
        }
        return p->a;  // Returns the reference to the matrix element.
    }
}

/**
 * Returns the element (i, j) of the matrix, or zero if the element does not exist.
 * This method is a getter: It does not change the matrix state.
 */
dcomplex SparseComplexMatrix::get(const int i, const int j) const {
    
    checkIndices(i, j);   // First check the validity of indices.
    checkSorted("get()"); // Then check that the matrix elements are sorted. Since get() is "const", it cannot call finalize itself.
    
    Triplet t;
    t.i = i; t.j = j;
    
    auto p = std::lower_bound(triplet.begin(), triplet.end(), t);  // Find the first element >= than the given triplet using binary search.
    
    if (p == triplet.end() || p->i != i || p->j != j) {// If the triplet does not exist in the current matrix.
        // Note that the first condition prevents the pointer "p" to be dereferenced out of range (because end() always points out of the vector).
        return dcomplex(0., 0.);
    }
    
    return p->a;
}

/**
 * Returns "true" is the matrix is symmetric, "false" otherwise.
 */
bool SparseComplexMatrix::isSymmetric() const {
    checkSorted("isSymmetric()");
    return symmetric;
}

/**
 * Returns "true" is the matrix is symmetric, "false" otherwise.
 * Argument "tol" is the tolerance over the relative difference.
 * Typically, it should be of the order of the machine epsilon (1e-16) or larger.
 * This method implicitly assumes that the matrix elements are already sorted and does not perform checking (this is why it should remain private).
 */
void SparseComplexMatrix::computeSymmetry() {
    
    if (nrow != ncol) {// If the matrix is non-square, then skip further calculations.
        symmetric = false;
        return;
    }
    const double tol = 1e-12;  // Tolerance over the relative error between A(i,j) and A(j,i).
    dcomplex aji;  // Symmetric element.
    
    //auto start_sym = std::chrono::steady_clock::now();
    
    symmetric = true; // The matrix is symmetric by default, until proven otherwise.
    
    for (const Triplet& t : triplet) {// Loop over the triplets.
        
        if (t.i > t.j) {// Only perform the check in the lower triangle.
            
            aji = get(t.j, t.i);  // Get the symmetric element.
            
            if (std::abs(t.a - aji) + std::abs(t.a - aji) > tol*(std::abs(t.a.real()) + std::abs(t.a.imag()) + 1.)) {// Symmetry criterion with tolerance.
                std::cout << TAG_WARN << "Element A(" << t.i << ", " << t.j << ") = " << t.a << " is different from A("
                          << t.j << ", " << t.i << ") = " << aji << "...\n";
                symmetric = false;
                return;
            }
        }
    }
    
    //double ctime_sym = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_sym).count();
    //std::cout << TAG_INFO << "computeSymmetry(): Time = " << ctime_sym << " s.\n";
}

/**
 * Sum the duplicate triplets (i.e., matrix elements at the same position (i,j)) of the sparse matrix after it is sorted.
 * Print a warning if duplicate entries are found.
 * This method implicitly assumes that the triplets are sorted in column-major ordering.
 */
void SparseComplexMatrix::sumDuplicate() {
    
    for (auto tp = triplet.end()-1; tp > triplet.begin(); tp--) {// Loop the vector backward (thus safely removing elements).
        if (*tp == *(tp-1)) {// If consecutive triplets have the same position (i,j) in the matrix.
            std::cout << TAG_WARN << "Summing triplets " << *tp << " and " << *(tp-1) << "...\n";
            (tp-1)->a += tp->a;
            triplet.erase(tp); // Remove the copied triplet.
            //std::cout << TAG_WARN << "Summed triplet is " << *(tp-1) << "...\n";
        }
    }
}

/**
 * Finalize the sparse matrix by sorting the matrix elements in column-major ordering.
 */
void SparseComplexMatrix::finalize() {
    if (not sorted) {// Only sort one time.
        std::sort(triplet.begin(), triplet.end());  // Sort is O(N*log(N)) in time.
        sumDuplicate(); // Sum the duplicate elements. About O(N) in time.
        sorted = true; // Declare the vector "triplet" as sorted before computing the symmetry.
        computeSymmetry(); // Determines if the matrix is symmetric, in principle O(N*log(N)) in time.
    }
}

/*******************************************************************************
 * PRINTING METHODS
 ******************************************************************************/

/**
 * Print essential information about the sparse matrix without showing the full content.
 */
void SparseComplexMatrix::printSummary(const std::string& name) const {
    std::cout << TAG_INFO << "Sparse complex matrix '" << name << "': Nrow=" << nrow << ", Ncol=" << ncol << ", Nnz=" << triplet.size() 
              << ", Capacity=" << triplet.capacity() << ", Density=" << 100.*density() << "%, Sorted=" << (sorted ? "true" : "false") 
              << ", Symmetric=" << (isSymmetric() ? "true" : "false") << ".\n";
}

/**
 * Prints the vector of triplets to std output.
 */
void SparseComplexMatrix::print(const std::string& name) const {
    std::cout << TAG_INFO << name << " (" << nrow << "x" << ncol << ", nnz=" << triplet.size() 
              << ", density=" << 100.*density() << "%, sorted=" << (sorted ? "true" : "false") << ") = [\n";
    for (const Triplet& t : triplet) {// Loop on the triplets.
        std::cout << "\t" << t << "\n";
    }
    std::cout << "]" << std::endl;
}

/**
 * Prints the sparsity pattern of the present sparse matrix to a PNG file.
 */
void SparseComplexMatrix::savePNG(const std::string& filename) const {
    
    Image img(nrow, ncol);
    img.fill(COLOR_WHITE); // Set white background.
    dcomplex elem;
    
    for (const Triplet& t : triplet) {// Loop on the nonzero elements.
        //img(t.i, t.j) = complexColor1(t.a);
        img(t.i, t.j) = complexColor2(t.a);
    }
    
    img.savePNG(filename);
}

/**
 * Prints the sparsity pattern of the present sparse matrix to a portable pixmap file, a PPM file (see: https://en.wikipedia.org/wiki/Netpbm).
 */
void SparseComplexMatrix::savePPM(const std::string& filename) const {
    
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
    
    for (const Triplet& t : triplet) {// Loop over the triplets.
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
    c.allocate(triplet.size()); // Preallocate the matrix to speedup the construction (avoid reallocation).
    
    for (const Triplet& t : triplet) {// Loop over the triplets.
        c(t.j, t.i) = std::conj(t.a);
    }
    c.finalize();  // Ensure that the returned matrix has proper ordering.
    
    return c;
}

/**
 * Multiplication of a sparse matrix by a scalar factor.
 */
SparseComplexMatrix operator*(const dcomplex scalar, const SparseComplexMatrix& a) {
    
    SparseComplexMatrix b(a.nrow, a.ncol);
    b.allocate(a.triplet.size()); // Preallocate the matrix to speedup the construction (avoid reallocation).
    
    for (const Triplet& t : a.triplet) {// Loop on the triplets of A.
        b(t.i, t.j) = scalar * t.a;
    }
    b.finalize(); // Ensure the matrix "B" is sorted (it should be so in principle).
    
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
    
    for (const Triplet& t : a.triplet) {// Loop over the triplets of A.
        
        for (int k = 0; k < b.getNcol(); k++) {// Loop over the columns of B.
            c(t.i, k) += t.a * b(t.j, k);   // C_ik = Sum_j A_ij B_jk
        }
        
    }
    return c;
}

/**
 * Temporary function to print arrays. Used only for testing purposes (otherwise deprecated).
 */
template <typename T>
void printArray(const std::string& name, const int n, const T* array) {
    std::cout << TAG_INFO << name << " = [";
    for (int i = 0; i < n; i++) {// Loop over the elements of the array.
        std::cout << " " << array[i] << " ";
    }
    std::cout << "]\n";
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
    a.checkSorted("solve*()");  // Check that the two sparse matrices are properly initialized (sorted).
    b.checkSorted("solve*()");
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
    
    //auto start_conversion = std::chrono::steady_clock::now();
    
    // 2. Convert the matrix into UMFPACK's compressed column format:
    int32_t nnz = a.triplet.size();  // Estimate the number of nonzero (with/without symmetry).
    auto pcol = new int32_t[a.ncol+1]();
    auto irow = new int32_t[nnz];
    auto real = new double[nnz];
    auto imag = new double[nnz];
    
    nnz = 0;  // Reset the number of nonzero to count them.
    int32_t jcol = 0;  // Current column index.
    
    for (const Triplet& t : a.triplet) {// Loop on the triplets in column-major order.
        
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
    
    //double ctime_conversion = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_conversion).count();
    //std::cout << TAG_INFO << "solveUmfpack(): Conversion time = " << ctime_conversion << " s.\n";
    //auto start_solution = std::chrono::steady_clock::now();
    
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
    umfpack_zi_numeric(pcol, irow, real, imag, symbolic, &numeric, control, info); // Perform numerical LU factorization.
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
    
    // Check for possible UMFPACK errors:
    if (info[UMFPACK_STATUS] != UMFPACK_OK) {
        switch (static_cast<int>(info[UMFPACK_STATUS])) {
            case UMFPACK_WARNING_singular_matrix:
                std::cout << TAG_WARN << "In solveUmfpack(): UMFPACK warning, singular matrix...\n";
                break;
            case UMFPACK_WARNING_determinant_underflow:
                std::cout << TAG_WARN << "In solveUmfpack(): UMFPACK warning, determinant underflow...\n";
                break;
            case UMFPACK_WARNING_determinant_overflow:
                std::cout << TAG_WARN << "In solveUmfpack(): UMFPACK warning, determinant overflow...\n";
                break;
            case UMFPACK_ERROR_invalid_matrix:
                std::cout << TAG_ERROR << "In solveUmfpack(): UMFPACK error, invalid matrix...\n";
                break;
            case UMFPACK_ERROR_invalid_Numeric_object:
                std::cout << TAG_ERROR << "In solveUmfpack(): UMFPACK error, invalid Numeric object, probably a memory overflow...\n";
                break;
            case UMFPACK_ERROR_invalid_Symbolic_object:
                std::cout << TAG_ERROR << "In solveUmfpack(): UMFPACK error, invalid Symbolic object, probably a memory overflow...\n";
                break;
            case UMFPACK_ERROR_out_of_memory:
                std::cout << TAG_ERROR << "In solveUmfpack(): UMFPACK error, not enough memory...\n";
                break;
            case UMFPACK_ERROR_different_pattern:
                std::cout << TAG_ERROR << "In solveUmfpack(): UMFPACK error. "
                    << "The pattern of the matrix has changed between the symbolic and numeric factorization...\n";
                break;
            case UMFPACK_ERROR_invalid_system:
                std::cout << TAG_ERROR << "In solveUmfpack(): UMFPACK error. "
                    << "The 'sys' argument provided to one of the solve routines is invalid...\n";
                break;
            case UMFPACK_ERROR_internal_error:
                std::cout << TAG_ERROR << "In solveUmfpack(): UMFPACK error, internal error (code=" << info[UMFPACK_STATUS] << ")...\n";
                break;
            default:
                std::cout << TAG_WARN << "In solveUmfpack(): UMFPACK returned an error code (code=" << info[UMFPACK_STATUS] << ")...\n";
        }
    }
    
    //double ctime_solution = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_solution).count();
    //std::cout << TAG_INFO << "solveUmfpack(): Solution time = " << ctime_solution << " s.\n";
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
    
    checkSolverInput(a, b, x);  // Check for possible invalid arguments.
    
    //====== Initialize a MUMPS instance ======
    ZMUMPS_STRUC_C id;  // Declare a new MUMPS object.
    id.comm_fortran = -987654;  // Tell MUMPS to use the default MPI communicator for Fortran (see MUMPS manual).
    id.par = 1;  // Declares that the host is involved in the factorization and solve phases (id.par=1 for sequential solve).
    id.sym = a.symmetric ? 2 : 0;  // Declares the symmetry of the matrix.
    id.job = -1; // Declares that this call initializes an instance of the package.
    zmumps_c(&id);
    
    //====== Enter the sparse matrix A ======
    #define ICNTL(I) icntl[(I)-1] // Macro such that indices match documentation (see MUMPS manual).
    id.ICNTL(5) = 0;  // Declare the use of assembled format (standard triplet format).
    const uint64_t a_nnz  = a.symmetric ? (a.triplet.size() + a.nrow)/2 : a.triplet.size();
    //std::cout << TAG_INFO << "In solveMumps(): Attempting to allocate a_nnz=" << a_nnz << std::endl;
    auto a_elem = new mumps_double_complex[a_nnz];
    auto a_irow = new int[a_nnz];
    auto a_jcol = new int[a_nnz];
    
    int nnz = 0;
    for (const Triplet& t : a.triplet) {// Loop on the triplets in column-major order.
        if (not a.symmetric || t.i >= t.j) {// If the matrix is symmetrix, then only retains the lower triangle (including diagonal).
            a_irow[nnz] = t.i + 1;  // Note the use of Fortran 1-based indices.
            a_jcol[nnz] = t.j + 1;
            a_elem[nnz] = mumps_double_complex{t.a.real(), t.a.imag()};
            nnz++;
        }
    }
    
    //std::cout << TAG_INFO << "solveMumps(): Allocated a_nnz=" << a_nnz << ", inserted nnz=" << nnz << "...\n";
    
    id.n   = a.nrow;
    id.nnz = a_nnz;
    id.irn = a_irow;
    id.jcn = a_jcol;
    id.a   = a_elem;
    
    //====== Enter the sparse right-hand side B ======
    const uint64_t b_nnz  = b.triplet.size();
    //std::cout << TAG_INFO << "In solveMumps(): Attempting to allocate b_nnz=" << b_nnz << std::endl;
    auto b_elem = new mumps_double_complex[b_nnz];
    auto b_irow = new int[b_nnz];
    auto b_pcol = new int[b.ncol+1];
    
    int jcol = 0; nnz = 0;
    b_pcol[0] = 1;
    for (const Triplet& t : b.triplet) {// Loop on the triplets in column-major order.
        if (t.j != jcol) {
            jcol = t.j;
            b_pcol[jcol] = nnz + 1;  // Note the use of Fortran 1-based indices.
        }
        b_irow[nnz] = t.i + 1;  // Note the use of Fortran 1-based indices.
        b_elem[nnz] = mumps_double_complex{t.a.real(), t.a.imag()};
        nnz++;
    }
    b_pcol[b.ncol] = nnz + 1;
    
    //printArray("b_pcol", b.ncol+1, b_pcol);
    
    id.nrhs        = b.ncol;
    id.lrhs        = b.nrow;
    id.nz_rhs      = b_nnz;
    id.rhs_sparse  = b_elem;
    id.irhs_sparse = b_irow;
    id.irhs_ptr    = b_pcol;
    
    //====== Solve the linear system A.X = B ======
    //std::cout << TAG_INFO << "In solveMumps(): Attempting to allocate x.getNrow() * x.getNcol() = " << uint64_t(x.getNrow() * x.getNcol()) << std::endl;
    id.rhs = new mumps_double_complex[uint64_t(x.getNrow() * x.getNcol())]; // Allocate space for the solution.
    //std::cout << TAG_INFO << "In solveMumps(): Allocation successful, now solving..." << std::endl;
    
    id.ICNTL(1)  = 0;  // Disable printing.
    id.ICNTL(2)  = 0;  // Disable printing.
    id.ICNTL(3)  = 0;  // Disable printing.
    id.ICNTL(4)  = 0;  // Disable printing.
    id.ICNTL(10) = 0;  // Declares no iterative refinement.
    id.ICNTL(16) = 1;  // Sets the number of OpenMP threads to 1.
    id.ICNTL(20) = 1;  // Declare sparse right-hand sides. Decision of exploiting sparsity of the RHS to accelerate the solution phase is done automatically.
    id.ICNTL(28) = 1;  // Declares the sequential computation of the ordering.
    
    id.job = 6; // Declares the computational phase (analysis, factorization, solve).
    zmumps_c(&id); // Call MUMPS solver (computationally intensive phase).
    
    if (id.infog[0] < 0) {// Check for possible errors.
        std::string msg = "In solveMumps(): MUMPS failure. INFOG(1)=" + std::to_string(id.infog[0])
                        + ", INFOG(2)=" + std::to_string(id.infog[1]) + ".";
        throw std::runtime_error(msg);
    }
    
    //====== Copy of the solution into X ======
    const int x_nrow = x.getNrow();
    mumps_double_complex elem;
    for (int j = 0; j < x.getNcol(); j++) {// Loop over the columns of X.
        for (int i = 0; i < x_nrow; i++) {// Loop over the rows of X.
            elem = id.rhs[i + j*x_nrow];
            x(i, j) = dcomplex(elem.r, elem.i);
        }
    }
    
    //====== Termination instructions ======
    delete[] id.rhs; // First ensure the memory allocated by the solution is free.
    id.job = -2; // Declares the termination of the MUMPS instance.
    zmumps_c(&id); // Termination.
    
    delete[] a_elem; // Finally, frees the memory.
    delete[] a_irow;
    delete[] a_jcol;
    delete[] b_elem;
    delete[] b_irow;
    delete[] b_pcol;
}

