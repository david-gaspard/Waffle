/****
 * @date Created on 2025-07-31 at 14:30:06 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the implementation of the WaveSystem methods.
 ***/
#include "WaveSystem.hpp"
#include "BaseTools.hpp"
#include <random>
#include <algorithm>
#include <iomanip>

/**
 * Constructor of the WaveSystem object.
 * 
 * Arguments:
 * 
 * sysname = String containing the name of the system (typically describing the system geometry) which is used to generate file output.
 * mesh    = Square mesh object.
 * kh      = Wavenumber times the lattice step, 2*pi*h/lambda. Also the phase accumulated across a lattice step (in radian).
 * density = Scatterer density per unit pixel. Between 0 and 1.
 * holscat = Lattice step divided by the scattering mean free path, h/lscat. Total length is not a well defined unit.
 * holabso = Lattice step divided by the absorption length, h/labso. Total length is not a well defined unit.
 */
WaveSystem::WaveSystem(const std::string& sysname, const SquareMesh& mesh, const double kh, const double density, const double holscat, const double holabso) : 
        sysname(sysname),  // Use initializer list to allocate matrix sizes.
        mesh(mesh),
        kh(checkWavenumber(kh)),
        density(checkDensity(density)),
        holscat(checkScattering(holscat)),
        holabso(checkAbsorption(holabso)),
        npoint(mesh.getNPoint()),
        ninput(mesh.getNBoundary(BND_INPUT)),
        noutput(mesh.getNBoundary(BND_OUTPUT)),
        ninputprop(computeNInputProp()),
        noutputprop(computeNOutputProp()),
        computed(false),
        hamiltonian(npoint, npoint),
        inputState(npoint, ninputprop),
        outputState(npoint, noutputprop),
        inputKlh(ninputprop, 1),
        outputKlh(noutputprop, 1),
        green(npoint, ninputprop)
    {
    //std::cout << TAG_INFO << "Creating WaveSystem...\n";
    computeHamiltonian();
    computeIOStates();
}

/***********************************************************
 * PRIVATE SETTERS
 ***********************************************************/

/**
 * Check the validity of the given wavenumber, k*h, and returns it.
 */
double WaveSystem::checkWavenumber(const double kh) const {
    if (kh < 0.) {
        throw std::invalid_argument("In checkWavenumber(): Wavenumber cannot be negative.");
    }
    else if (kh > 2.) {
        throw std::invalid_argument("In checkWavenumber(): Wavenumber cannot exceed Nyquist-Shannon bound (kh<2).");
    }
    return kh;
}

/**
 * Check the validity of the given density, and returns it.
 */
double WaveSystem::checkDensity(const double density) const {
    if (density < 0. || density > 1.) {
        throw std::invalid_argument("In checkDensity(): Density must be between 0 and 1.");
    }
    return density;
}

/**
 * Assigns the scattering strength, h/lscat.
 */
double WaveSystem::checkScattering(const double holscat) const {
    if (holscat < 0.) {
        throw std::invalid_argument("In checkScattering(): Scattering strength cannot be negative.");
    }
    else if (kh < holscat) {
        std::cout << TAG_WARN << "Scattering strength is large (k*lscat=" << kh/holscat << " <1). Localization may occur.\n";
    }
    return holscat;
}

/**
 * Assigns the absorption strength, h/labso.
 */
double WaveSystem::checkAbsorption(const double holabso) const {
    return holabso;
}

/***********************************************************
 * PUBLIC GETTERS
 ***********************************************************/

/**
 * Returns the total number of points in the mesh.
 */
int WaveSystem::getNPoint() const {
    return npoint;
}

/**
 * Returns the i-th point of the mesh.
 */
MeshPoint WaveSystem::getPoint(const unsigned int ipoint) const {
    return mesh.getPoint(ipoint);
}

/**
 * Returns the number of input points, also the number of input channels.
 */
int WaveSystem::getNInput() const {
    return ninput;
}

/**
 * Returns the number of output points, also the number of output channels.
 */
int WaveSystem::getNOutput() const {
    return noutput;
}

/**
 * Returns the number of propagating modes in the input lead(s).
 */
int WaveSystem::getNInputProp() const {
    return ninputprop;
}

/**
 * Returns the number of propagating modes in the output lead(s).
 */
int WaveSystem::getNOutputProp() const {
    return noutputprop;
}

/**
 * Returns the wavenumber times the lattice step, k*h = 2*pi*h/lambda.
 */
double WaveSystem::getWavenumber() const {
    return kh;
}

/**
 * Returns the density of scatterers per unit pixel.
 */
double WaveSystem::getDensity() const {
    return density;
}

/**
 * Returns the value of "holscat" which is defined by h/lscat, where "h" is the lattice step
 * (the unit length) and "lscat" is the scattering mean free path.
 */
double WaveSystem::getScattering() const {
    return holscat;
}

/**
 * Returns the value of "holabso" which is defined by h/labso, where "h" is the lattice step
 * (the unit length) and "labso" is the ballistic absorption length.
 */
double WaveSystem::getAbsorption() const {
    return holabso;
}

/**
 * Returns a deep copy of the name of the UsadelSystem.
 */
std::string WaveSystem::getName() const {
    return sysname;
}

/***********************************************************
 * PRINTING METHODS
 ***********************************************************/

/**
 * Returns a unique output filename for the given "dataname" with given file "extension" (with dot).
 */
std::string WaveSystem::uniqueFile(const std::string& dataname, const std::string& extension) const {
    std::string filename;
    uniqueFilename("out/" + sysname + "/" + dataname + "/" + dataname + "_", extension, filename);
    return filename;
}

/**
 * Returns the summary of essential information about the current wave system.
 */
std::vector<std::string> WaveSystem::summary() const {
    std::vector<std::string> smr;
    smr.push_back("WaveSystem with sysname='" + sysname + "', Npoint=" + std::to_string(npoint) + ", kh=" + std::to_string(kh)
        + ", lambda/h=[" + std::to_string(PI/std::asin(kh/2)) + ", " + std::to_string(PI/(std::sqrt(2.)*std::asin(kh/std::sqrt(8.)))) + "]");
    smr.push_back("density=" + std::to_string(density) + ", h/lscat=" + std::to_string(holscat) + ", h/labso=" + std::to_string(holabso)
        + ", Ninputprop/Ninput=" + std::to_string(ninputprop) + "/" + std::to_string(ninput) 
        + ", Noutputprop/Noutput=" + std::to_string(noutputprop) + "/" + std::to_string(noutput));
    smr.push_back("DOSinput=" + std::to_string(dosinput) + ", DOSoutput=" + std::to_string(dosoutput) + ", DOSlattice=" + std::to_string(doslattice)
        + ", DOSfree=" + std::to_string(DOSFREE) + ", Hamiltonian_sparse_density=" + std::to_string(100.*hamiltonian.density()) + "%");
    return smr;
}

/**
 * Print essential information about the current wave system to standard output.
 */
void WaveSystem::printSummary() const {
    for (const std::string& line : summary()) {// Loop over the lines of the summary.
        std::cout << TAG_INFO << line << "\n";
    }
}

/**
 * Prints essential information about the sparse Hamiltonian matrix.
 */
void WaveSystem::infoHamiltonian() const {
    hamiltonian.printSummary("Hamiltonian");
}

/**
 * Prints the sparsity pattern of the Hamiltonian to a portable pixmap file, a PPM file (see: https://en.wikipedia.org/wiki/Netpbm).
 */
void WaveSystem::plotHamiltonian() const {
    
    const std::string filename = uniqueFile("hamiltonian", ".png");
    const auto start = std::chrono::steady_clock::now();
    
    hamiltonian.savePNG(filename);
    
    double ctime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start).count();
    std::cout << TAG_INFO << "Hamiltonian saved to file '" << filename << "' in " << ctime << " s.\n";
}

/**
 * Prints the sparsity pattern of the input state matrix to a portable pixmap file, a PPM file (see: https://en.wikipedia.org/wiki/Netpbm).
 */
void WaveSystem::plotInputState() const {
    
    const std::string filename = uniqueFile("input_state", ".png");
    const auto start = std::chrono::steady_clock::now();
    
    inputState.savePNG(filename);
    
    double ctime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start).count();
    std::cout << TAG_INFO << "Input state saved to file '" << filename << "' in " << ctime << " s.\n";
}

/**
 * Prints the sparsity pattern of the output state matrix to a portable pixmap file, a PPM file (see: https://en.wikipedia.org/wiki/Netpbm).
 */
void WaveSystem::plotOutputState() const {
    
    const std::string filename = uniqueFile("output_state", ".png");
    const auto start = std::chrono::steady_clock::now();
    
    outputState.savePNG(filename);
    
    double ctime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start).count();
    std::cout << TAG_INFO << "Input state saved to file '" << filename << "' in " << ctime << " s.\n";
}

/**
 * Plot the mesh contained in the present wave system.
 */
void WaveSystem::plotMesh() const {
    mesh.plotMesh(uniqueFile("mesh", ".csv"));
}

/**
 * Save the given intensity profile to a CSV file and plot it using an external script.
 */
void WaveSystem::plotIntensity(const RealMatrix& intensity, const std::string& description, const std::string& dataname) const {
    
    if (intensity.getNrow() != npoint) {// First check for possible errors:
        std::string msg = "In plotIntensity(): Invalid intensity matrix, received nrow=" + std::to_string(intensity.getNrow()) 
                        + ", expected nrow=" + std::to_string(npoint) + ".";
        throw std::invalid_argument(msg);
    }
    else if (description.empty()) {
        throw std::invalid_argument("In plotIntensity(): Empty 'description'. Please provide a description of the intensity profile.");
    }
    else if (dataname.empty()) {
        throw std::invalid_argument("In plotIntensity(): Empty 'dataname'. Please provide a short folder name for the intensity profile.");
    }
    
    const int nstate = intensity.getNcol();  // Number of states in the wavefunction.
    const char* sep = ", ";  // Separator used between entries of the CSV file.
    const int prec = 16;     // Precision used in printing double precision values.
    const std::string filename = uniqueFile(dataname, ".csv");
    const double fsize = ( (prec+4.)*nstate + 3.*DIMENSION*(std::log10(npoint)+2.) ) * static_cast<double>(npoint);  // Roughly estimated file size in bytes (octets).
    std::cout << TAG_INFO << "Save intensity to file '" << filename << "', size ~" << (fsize/1e6) << " Mo...\n";
    
    std::ofstream ofs(filename); // Open the output file.
    ofs << std::setprecision(prec); // Set the printing precision.
    writeTimestamp(ofs, "%% "); // Apply a timestamp at the beginning.
    
    for (const std::string& line : summary()) {// Write the summary to the file header.
        ofs << "%% " << line << "\n";
    }
    ofs << "%% dataname='" + dataname + "'\n%% Info: " << description << "\n"
        << "x" << sep << "y" << sep << "north" << sep << "south" << sep << "east" << sep << "west";
    
    for (int istate = 0; istate < nstate; istate++) {// Loop over the input modes to finish the column names.
        ofs << sep << "I" << istate;
    }
    ofs << "\n";
    
    MeshPoint p;
    dcomplex psi;
    
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
        
        p = mesh.getPoint(ipoint); // Extract the point to get its coordinates.
        ofs << p.x << sep << p.y << sep
            << boundaryTypeString(p.north) << sep << boundaryTypeString(p.south) << sep 
            << boundaryTypeString(p.east)  << sep << boundaryTypeString(p.west);
        
        for (int istate = 0; istate < nstate; istate++) {// Loop over the input modes.
            ofs << sep << intensity(ipoint, istate); // Save the square modulus to the file.
        }
        ofs << "\n";
    }
    ofs.close();  // Close the stream before calling an external script (this may cause I/O trouble).
    
    const std::string cmd = "plot/plot_map.py lin I0 " + std::to_string(holscat) +  " " + filename;
    std::cout << TAG_EXEC << cmd << "\n";
    if (std::system(cmd.c_str())) {
        std::cout << TAG_WARN << "The plot script returned an error.\n";
    }
}

/**
 * Save the square modulus of the Green functions associated to each input mode to a CSV file and
 * call an external script to plot the lowest mode.
 */
void WaveSystem::plotGreenFunction() {
    
    computeGreenFunction();  // Ensure that the Green function has already been computed.
    
    // Compute the square modulus of the Green function:
    RealMatrix intensity(npoint, ninputprop);
    dcomplex psi;
    for (int imode = 0; imode < ninputprop; imode++) {// Loop over the modes.
        for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
            psi = green(ipoint, imode);
            intensity(ipoint, imode) = psi.real()*psi.real() + psi.imag()*psi.imag();
        }
    }
    
    std::string description = "Modal Green functions";
    plotIntensity(intensity, description, "green");
}

/**
 * Save the square modulus of the first transmission eigenstates to a CSV file and call an external script to plot them.
 * "nstate" is the number of desired transmission eigenstates, starting from the largest transmission value.
 */
void WaveSystem::plotTransmissionStates(const int nstate) {
    
    // 1. First check for possible invalid values of "nstate":
    if (nstate <= 0 || nstate > std::min(noutputprop, ninputprop)) {
        std::string msg = "In plotTransmissionStates(): Invalid number of transmission eigenstates, received nstate=" + std::to_string(nstate) 
                        + ", expected in 1.." + std::to_string(std::min(noutputprop, ninputprop)) + ".";
        throw std::invalid_argument(msg);
    }
    
    // 2. Compute the transmission eigenstates:
    ComplexMatrix tstate(npoint, nstate);
    RealMatrix tval(nstate, 1);
    computeTransmissionStates(tstate, tval);
    
    // 2. Then compute the square modulus of the transmission eigenstates:
    RealMatrix intensity(npoint, nstate);
    dcomplex psi;
    for (int istate = 0; istate < nstate; istate++) {// Loop over the transmission eigenstates.
        for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
            psi = tstate(ipoint, istate);
            intensity(ipoint, istate) = psi.real()*psi.real() + psi.imag()*psi.imag();
        }
    }
    
    std::string description = "Transmission eigenstates with Tval=[";
    for (int istate = 0; istate < nstate; istate++) {// Loop over the transmission eigenvalues.
        description += std::to_string(tval(istate, 0)) + ", ";
    }
    description += "] (same order).";
    plotIntensity(intensity, description, "tstate");
}

/***********************************************************
 * COMPUTATIONAL METHODS
 ***********************************************************/

/**
 * Creates the free part of the "Hamiltonian", (d_x^2 + d_y^2 + k^2) * h^2, and store the result within the present WaveSystem object.
 */
void WaveSystem::computeHamiltonian() {
    std::cout << TAG_INFO << "Building the Hamiltonian... ";
    const auto start_build = std::chrono::steady_clock::now(); // Gets the current time.
    
    const dcomplex kh2 = dcomplex(kh*kh, MEPS);
    const std::vector<Opening> opening = mesh.getOpening(); // Extract the list of openings.
    Opening op;
    MeshPoint p;
    uint iop;
    int np, i, j, ip, jp;
    
    // 1. Preallocate the sparse Hamiltonian and precompute the opening matrices:
    int nnzub = 5*npoint; // Computes an upper bound on the number of nonzero elements.
    std::vector<ComplexMatrix> opmatrix; // Precompute the opneing matrices.
    for (const Opening& op : opening) {// Loop over the openings.
        np = op.index.size();
        nnzub += (np - 1)*np;
        opmatrix.push_back(openingMatrix(kh2, np));
    }
    hamiltonian.allocate(nnzub);
    //std::cout << TAG_INFO << "Preallocated Hamiltonian with nnzub = " << nnzub << "\n";
    
    // 2. Build the Hamiltonian row per row:
    for (i = 0; i < npoint; i++) {//Loop over the points of the mesh.
        
        p = mesh.getPoint(i);
        
        if (not p.isOpening()) {// If the point is not in an opening.
            
            hamiltonian(i, i) = -4. + kh2;  // Add the diagonal element of the Hamiltonian H(i, i) = -4 + (k*h)^2.
            
            for (const Direction dir : allDirections) {// Loop over all directions.
                j = p.neighbor(dir);  // Extract the index of the neighboring point.
                if (j >= 0) {// If the neigbhoring point is in the mesh, then add H(i, j) = 1.
                    hamiltonian(i, j) = 1.;
                }
            }
        }
        else {// If the point is actually in an opening.
            // Find to which opening the point p(i) belongs:
            for (iop = 0; iop < opening.size(); iop++) {// Loop over the openings.
                op = opening.at(iop);
                auto ptr = lower_bound(op.index.begin(), op.index.end(), i);  // Find the point p(i) in the opening "op".
                if (ptr != op.index.end() && *ptr == i) {// If the present opening contains the point p(i).
                    ip = std::distance(op.index.begin(), ptr); // Get the index of p(i) inside the opening.
                    np = op.index.size(); // Number of points in this opening.
                    for (jp = 0; jp < np; jp++) {// Loop over the points in the present opening.
                        hamiltonian(i, op.index.at(jp)) = opmatrix.at(iop)(ip, jp);
                    }
                }
            }
            for (const Direction dir : allDirections) {// Loop over the directions.
                j = p.neighbor(dir);
                if (j >= 0 && not mesh.getPoint(j).isOpening()) {// If point of index "j" belongs to the mesh, and is not an opening, then add H(i, j) = 1.
                    hamiltonian(i, j) = 1.;
                }
            }
        }
    }
    
    // 3. Finalize the sparse Hamiltonian (sort the matrix elements in column-major ordering):
    hamiltonian.finalize();
    
    // Print some warnings:
    if (hamiltonian.getNnz() > nnzub) {
        std::cout << TAG_WARN << "Insufficient matrix preallocation. Added nnz=" << hamiltonian.getNnz() 
                  << " elements, but allocated only nnzub=" << nnzub << ".\n";
    }
    
    // Measure the build time for information:
    double ctime_build = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_build).count();
    std::cout << "Done in " << ctime_build << " s.\n";
}

/**
 * Assigns a realization of the disorder to the Hamiltonian.
 */
void WaveSystem::setDisorder(const uint64_t seed) {
    //std::cout << TAG_INFO << "Setting disorder with seed=" << seed << "..." << std::endl;
    
    computed = false;  // Each time the Hamiltonian is modified, the Green function is no more valid and must be computed again.
    double uh2;
    dcomplex kh2 = dcomplex(kh*kh, MEPS);
    
    // Standard deviation of the random potential uh2 = U(x, y) * h^2 approximately producing the scattering strength h/lscat : 
    const double sigma = std::sqrt(kh*holscat/(PI*doslattice*density));
    //const double sigma = std::sqrt(kh*holscat/(PI*doslattice));
    
    std::mt19937_64 rng;  // Instantiate the standard Mersenne Twister random number generator (64-bit-return version).
    rng.seed(seed);       // Initialize the random generator with the given seed.
    std::normal_distribution<double> random_normal(0., sigma);
    std::uniform_real_distribution<double> random_uniform(0., 1.);
    
    for (int i = 0; i < npoint; i++) {//Loop over the points of the mesh.
        
        if (not mesh.getPoint(i).isOpening()) {// If the point is not in an opening, then modifies the potential.
            if (random_uniform(rng) < density) {// If the pixel is a scattering pixel (probability given by 'density').
                uh2 = random_normal(rng);  // Generate a Gaussian random number with standard deviation "sigma".
            }
            else {// Otherwise, then reset the value of the potential to zero.
                uh2 = 0.;
            }
            hamiltonian(i, i) = -4. + kh2 - uh2;  // In any case, replace the diagonal element of the Hamiltonian H(i, i) = -4 + (k*h)^2 - U*h^2.
        }
    }
}

/**
 * Computes and return the number of input propagating modes.
 */
int WaveSystem::computeNInputProp() const {
    
    const double fac = (2./PI) * std::asin(kh/2. - KLHMIN*KLHMIN/(4.*kh));
    int np, ninputprop = 0; // Initialize the number of input propagating modes (it will be incremented).
    
    for (const Opening& op : mesh.getOpening()) {// Loop over the openings.
        
        if (op.bndtype == BND_INPUT) {// If the opening is an input, then construct the modes.
            
            np = op.index.size();  // Number of point in the current input.
            
            /**
             * kh2 + d2ev > KLHMIN^2 (> 0)
             * 
             * d2ev = -4.*sin((i+1)*PI/(2*(n+1)))^2
             * 
             * kh2 - 4*sin((i+1)*PI/(2*(n+1)))^2 > KLHMIN^2, i=[0, n-1]
             * 
             * ninputprop < (2*(n+1)/PI) * asin(sqrt((kh2 - KLHMIN^2)/4))
             * 
             * Using the approx: sqrt((kh2 - KLHMIN^2)/4) <= kh/2 - KLHMIN^2/(4*kh)
             * 
             * ninputprop = std::floor( (2*(n+1)/PI) * std::asin(kh/2 - KLHMIN^2/(4*kh)) );
             */
            
            ninputprop += std::floor(fac*(np + 1));
        }
    }
    
    if (ninputprop == 0) {// Check for possible errors.
        throw std::invalid_argument("In computeNInputProp(): The number of input progpating modes is zero. Please increase the number of points in the opening.");
    }
    
    return ninputprop;
}

/**
 * Computes and return the number of input propagating modes.
 */
int WaveSystem::computeNOutputProp() const {
    
    const double fac = (2./PI) * std::asin(kh/2. - KLHMIN*KLHMIN/(4.*kh));
    int np, noutputprop = 0; // Initialize the number of input propagating modes (it will be incremented).
    
    for (const Opening& op : mesh.getOpening()) {// Loop over the openings.
        
        if (op.bndtype == BND_OUTPUT) {// If the opening is an input, then construct the modes.
            np = op.index.size();  // Number of point in the current input.
            noutputprop += std::floor(fac*(np + 1));
        }
    }
    
    if (noutputprop == 0) {// Check for possible errors.
        throw std::invalid_argument("In computeNInputProp(): The number of input progpating modes is zero. Please increase the number of points in the opening.");
    }
    
    return noutputprop;
}

/**
 * Compute the input and output matrices containing the input and output modes.
 */
void WaveSystem::computeIOStates() {
    std::cout << TAG_INFO << "Building the input/output states... ";
    const auto start_build = std::chrono::steady_clock::now(); // Gets the current time.
    
    // Construct the input/output modes :
    const dcomplex kh2 = dcomplex(kh*kh, MEPS);
    dcomplex d2pkh2, klh;
    const double klhmin2 = KLHMIN*KLHMIN; // Minimum value of klh^2 for a mode to be considered as a propagating.
    int np;          // Number of points in the current opening.
    int jinput = 0;  // Current number of input modes.
    int joutput = 0; // Current number of output modes.
    dosinput = 0.;   // Initialize the density of states in the input lead(s).
    dosoutput = 0.;  // Initialize the density of states in the output lead(s).
    
    for (const Opening& op : mesh.getOpening()) {// Loop over the openings.
        
        if (op.bndtype == BND_INPUT) {// If the opening is an input, then construct the modes.
            
            np = op.index.size();  // Number of point in the current input.
            ComplexMatrix u = modalMatrix(np);   // Construct the matrix of modes (each mode is normalized to 1).
            
            for (int j = 0; j < np; j++) {// Loop over the modes of U (columns of U).
                
                d2pkh2 = laplacianEigenvalue(j, np) + kh2;  // j^th eigenvalue of D_y^2 + (kh)^2 for an operator of size "np". Close to klh^2.
                
                if (d2pkh2.real() > klhmin2) {// If the mode is propagating.
                    
                    klh = std::sqrt( d2pkh2 * ( 1. - d2pkh2/4. ) );  // Compute the effective longitudinal wavenumber.
                    
                    for (int i = 0; i < np; i++) {// Loop over the points in the opening (rows of U).
                        inputState(op.index.at(i), jinput) = u(i, j);
                    }
                    inputKlh(jinput, 0) = klh;
                    dosinput += 1./klh.real();  // Increments the DOS in the input lead(s).
                    jinput++;
                }
            }
        }
        else if (op.bndtype == BND_OUTPUT) {// If the opening is an output, then construct the modes.
            
            np = op.index.size();  // Number of point in the current input.
            ComplexMatrix u = modalMatrix(np);   // Construct the matrix of modes (each mode is normalized to 1).
            
            for (int j = 0; j < np; j++) {// Loop over the modes of U (columns of U).
                
                d2pkh2 = laplacianEigenvalue(j, np) + kh2;  // j^th eigenvalue of D_y^2 + (kh)^2 for an operator of size "np". Close to klh^2.
                
                if (d2pkh2.real() > klhmin2) {// If the mode is propagating.
                    
                    klh = std::sqrt( d2pkh2 * ( 1. - d2pkh2/4. ) );  // Compute the effective longitudinal wavenumber.
                    
                    for (int i = 0; i < np; i++) {// Loop over the points in the opening (rows of U).
                        outputState(op.index.at(i), joutput) = u(i, j);
                    }
                    outputKlh(joutput, 0) = klh;
                    dosoutput += 1./klh.real();  // Increments the DOS in the output lead(s).
                    joutput++;
                }
            }
        }
    }
    
    // Finalize the sparse matrices (sort the elements):
    inputState.finalize();
    outputState.finalize();
    
    // Print some warnings:
    if (jinput != ninputprop) {
        std::cout << TAG_WARN << "Computed jinput=" << jinput << " input states, but planned ninputprop=" << ninputprop << ".\n";
    }
    if (joutput != noutputprop) {
        std::cout << TAG_WARN << "Computed joutput=" << joutput << " output states, but planned noutputprop=" << noutputprop << ".\n";
    }
    
    // Normalize the density of states:
    dosinput  /= 2*PI*ninput;  // Note that it is the total number of input/output modes (including evanescent modes).
    dosoutput /= 2*PI*noutput;
    
    // Compute the exact free density of states on a square lattice:
    const double mkh = 1. - kh*kh/4.;
    doslattice = ellipticK(1. - 1./(mkh*mkh))/(2.*PI*PI*mkh);
    
    if (std::abs(dosinput - doslattice) > 0.05*doslattice) {// Tolerance of 5% on the DOS.
        std::cout << TAG_WARN << "DOSinput=" << dosinput << " is different from DOSlattice=" << doslattice << ", meaning that an input lead resonates (inputKlh.real.min=" << inputKlh.real().min() << "). You may consider changing the wavenumber...\n";
    }
    if (std::abs(dosoutput - doslattice) > 0.05*doslattice) {// Tolerance of 5% on the DOS.
        std::cout << TAG_WARN << "DOSoutput=" << dosoutput << " is different from DOSlattice=" << doslattice << ", meaning that an output lead resonates (outputKlh.real.min=" << outputKlh.real().min() << "). You may consider changing the wavenumber...\n";
    }
    
    // Measure the build time for information:
    double ctime_build = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_build).count();
    std::cout << "Done in " << ctime_build << " s.\n";
}

/**
 * Compute the retarded Green function between a point in the input lead(s) to a point in the output lead(s).
 * If the Green function has already been computed and the Hamiltonian has not changed, this method does nothing.
 * This method contains the most time-consuming operation in the program.
 */
void WaveSystem::computeGreenFunction() {
    if (not computed) {// This function only does the computation if the flag "computed" is "false".
        //std::cout << TAG_INFO << "Solving the sparse system now..." << std::endl;
        solveUmfpack(hamiltonian, inputState, green);
            // Solve the system using UMFPACK (no parallelization, no iterative refinement).
        //solveMumps(hamiltonian, inputState, green);
            // Solve the system using MUMPS sequentially (MPI disabled, no iterative refinement).
            // Note that MUMPS is faster than UMFPACK but uses more memory.
        computed = true;  // Declare the Green function as computed.
    }
}

/**
 * Compute the transmission matrix "tmat" from the Green function using the Fisher & Lee relation.
 * t_ij = 2*I*sqrt(k_i k_j) G(Output mode i | Input mode j). Size: (noutput, ninput).
 */
void WaveSystem::transmissionMatrix(ComplexMatrix& tmat) {
    
    if (tmat.getNrow() != noutputprop || tmat.getNcol() != ninputprop) {// First check for possible errors.
        std::string msg = "In transmissionMatrix(): Invalid transmission matrix size, received (" 
                        + std::to_string(tmat.getNrow()) + ", " + std::to_string(tmat.getNcol()) + "), expected (" 
                        + std::to_string(noutputprop) + ", " + std::to_string(ninputprop) + ").";
        throw std::invalid_argument(msg);
    }
    
    computeGreenFunction();  // Ensure that the Green function has been computed (this does nothing if it is so).
    tmat = outputState.conj() * green;  // Project the Green function over the output states.
    
    // Apply the Fisher & Lee relation :
    for (int i = 0; i < noutputprop; i++) {// Loop over the output channels (rows).
        for (int j = 0; j < ninputprop; j++) {// Loop over the input channels (columns).
            tmat(i, j) = 2. * I * std::sqrt(inputKlh(j, 0).real() * outputKlh(i, 0).real()) * tmat(i, j);
            // Note that the real parts strip the evanescent modes.
        }
    }
}

/**
 * Compute the reflection matrix "rmat" from the Green function using hte Fisher & Lee relation.
 * r_ij = -delta_ij + 2*I*sqrt(k_i k_j) G(Input mode i | Input mode j). Size: (ninput, ninput).
 */
void WaveSystem::reflectionMatrix(ComplexMatrix& rmat) {
    
    if (rmat.getNrow() != ninputprop || rmat.getNcol() != ninputprop) {// First check for possible errors.
        std::string msg = "In reflectionMatrix(): Invalid reflection matrix size, received (" 
                        + std::to_string(rmat.getNrow()) + ", " + std::to_string(rmat.getNcol()) + "), expected (" 
                        + std::to_string(ninputprop) + ", " + std::to_string(ninputprop) + ").";
        throw std::invalid_argument(msg);
    }
    
    computeGreenFunction();  // Ensure that the Green function has been computed (this does nothing if it is so).
    rmat = inputState.conj() * green;  // Project the Green function over the input states.
    
    // Apply the Fisher & Lee relation :
    for (int i = 0; i < ninputprop; i++) {// Loop over the input channels (rows).
        for (int j = 0; j < ninputprop; j++) {// Loop over the input channels (columns).
            rmat(i, j) = 2. * I * std::sqrt(inputKlh(j, 0).real() * inputKlh(i, 0).real()) * rmat(i, j);
            // Note that the real parts strip the evanescent modes.
        }
        rmat(i, i) -= 1.;  // Note Fisher & Lee: r_ij = -delta_ij + 2*I*sqrt(k_i*k_j)*G(Input mode i | Input mode j).
    }
}

/**
 * Check that the residual of the solution of the linear system is reasonably close to zero.
 */
void WaveSystem::checkResidual() {
    
    computeGreenFunction();  // Ensure that the Green function has been computed (this does nothing if it is so).
    
    std::cout << TAG_INFO << "Computing the residual...\n";
    const ComplexMatrix inputState_product = hamiltonian * green; // Recompute the input state from the solution of the linear system.
    
    // Compare the matrices elementwise:
    const double tol = 1e-11; // Tolerance over the relative error (elementwise).
    double res, res_total = 0.;
    dcomplex elem, elem_expc;
    
    for (int i = 0; i < npoint; i++) {// Loop over the mesh points.
        for (int j = 0; j < ninputprop; j++) {// Loop over the input modes.
            elem = inputState_product(i, j);
            elem_expc = inputState.get(i, j);
            res = std::abs(elem - elem_expc);
            res_total += res;
            if (res > tol*(std::abs(elem_expc) + 1.)) {// If the residual is two large.
                std::cout << TAG_WARN << "Matrix element (" << i << ", " << j << ") is "
                          << elem << ", expected " << elem_expc << " (diff=" << res << ").\n";
            }
        }
    }
    std::cout << TAG_INFO << "Average residual = " << (res_total/npoint)/ninputprop << " (total=" << res_total << ").\n";
}

/**
 * Check the unitary of the propagation, i.e., check that in the absence of absorption, the relation associated to probability conservation,
 * t*t.conj() + r*r.conj() = identityMatrix(), is verified. This method does the checking for several realizations of the disorder (number "nseed").
 * This method is mainly used for testing purposes.
 * If "showtval" is "true", then compute the transmission eigenvalues and print them to standard output.
 */
void WaveSystem::checkUnitarity(const bool showtval) {
    ComplexMatrix tmat(noutputprop, ninputprop), rmat(ninputprop, ninputprop), unitarity(ninputprop, ninputprop);
    transmissionMatrix(tmat);
    reflectionMatrix(rmat);
    
    unitarity = tmat.conj() * tmat + rmat.conj() * rmat - identityMatrix(ninputprop);  // The unitarity matrix should ideally be zero on output.
    
    std::cout << TAG_INFO << "Unitarity: t^H * t + r^H * r - 1 = " << unitarity.norm() << "\n";
    
    if (showtval) {// Compute and print the transmission eigenvalues.
        const int ntval = std::min(noutputprop, ninputprop);  // Number of transmission eigenvalues.
        ComplexMatrix u(noutputprop, noutputprop), vh(ninputprop, ninputprop);
        RealMatrix tval(ntval, 1);
        svd(tmat, tval, u, vh);  // Compute the singular value decomposition (SVD).
        for (int i = 0; i < ntval; i++) tval(i, 0) = tval(i, 0)*tval(i, 0);  // COnvert singular values of "t" to transmission eigenvalues.
        tval.transpose().print("Tval");
    }
}

/**
 * Compute the intensity corresponding to a isotropic mixture of all input modes, and add it to the given data matrix.
 * The intensity is normalized to get unit incident intensity.
 */
void WaveSystem::addIIsotropic(RealMatrix& iavg) {
    
    // 1. First check for possible errors:
    if (iavg.getNrow() != npoint || iavg.getNcol() != 1) {
        std::string msg = "In addIIsotropic(): Invalid size of 'iavg'. Received ("
                        + std::to_string(iavg.getNrow()) + ", " + std::to_string(iavg.getNcol()) + "), expected ("
                        + std::to_string(npoint) + ", 1).";
        throw std::invalid_argument(msg);
    }
    
    computeGreenFunction(); // Ensure that the Green function is computed.
    
    // Add the mode-averaged intensity:
    double intensity;
    dcomplex psi;
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
        
        // Compute the mode-averaged intensity at the current position:
        intensity = 0.;
        for (int imode = 0; imode < ninputprop; imode++) {// Loop over the input modes to construct the transmission eigenstate at "ipoint".
            psi = green(ipoint, imode) * I*std::sqrt(2.*ninputprop/(PI*dosinput)) * std::sqrt(inputKlh(imode, 0).real());
            intensity += psi.real()*psi.real() + psi.imag()*psi.imag();
        }
        intensity /= ninputprop;
        iavg(ipoint, 0) += intensity;
        // Warning: this formula does not take into account multiple input leads.
    }
}

/**
 * Compute the intensity corresponding to the given input mode "imode", and add it to the given data matrix.
 * The intensity is normalized to get unit incident intensity.
 */
void WaveSystem::addIMode(RealMatrix& iavg, const int imode) {
    
    // 1. First check for possible errors:
    if (iavg.getNrow() != npoint || iavg.getNcol() != 1) {
        std::string msg = "In addIMode(): Invalid size of 'iavg'. Received ("
                        + std::to_string(iavg.getNrow()) + ", " + std::to_string(iavg.getNcol()) + "), expected ("
                        + std::to_string(npoint) + ", 1).";
        throw std::invalid_argument(msg);
    }
    else if (imode < 0 || imode >= ninputprop) {
        std::string msg = "In addIMode(): Invalid mode index 'imode'. Received "
                        + std::to_string(imode) + ", expected in 0.." + std::to_string(ninputprop-1) + ".";
        throw std::invalid_argument(msg);
    }
    
    computeGreenFunction(); // Ensure that the Green function is computed.
    
    // Add the mode-averaged intensity:
    dcomplex psi;
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
        psi = green(ipoint, imode) * 2.*I*std::sqrt(ninput) * inputKlh(imode, 0).real();
        iavg(ipoint, 0) += psi.real()*psi.real() + psi.imag()*psi.imag();
    }
    // Warning: this formula does not take into account multiple input leads.
}

/**
 * Compute all the transmission eigenvalues and their associated transmission eigenstates.
 * 
 * Arguments:
 * 
 * tstate  = On output, selected transmission eigenstates. Size: (npoint, nstate).
 * tval    = On output, list of transmission eigenvalues. Size: (nstate, 1).
 * Note that "nstate" is the desired number of transmission eigenstates, starting from the largest transmission eigenvalue.
 */
void WaveSystem::computeTransmissionStates(ComplexMatrix& tstate, RealMatrix& tval) {
    
    // 1. First check for possible errors:
    const int nstate = tval.getNrow(); // Desired number of transmission eigenstates (must be < nstate_max).
    const int nstate_max = std::min(ninputprop, noutputprop); // Expected maximum number of transmission eigenstates.
    
    if (nstate <= 0 || nstate > nstate_max) {
        std::string msg = "In computeTransmissionStates(): Invalid number of transmission eigenstates. Received "
                        + std::to_string(nstate) + ", expected in 1.." + std::to_string(nstate_max) + ".";
        throw std::invalid_argument(msg);
    }
    else if (tstate.getNrow() != npoint || tstate.getNcol() != nstate) {
        std::string msg = "In computeTransmissionStates(): Invalid size of 'tstate'. Received ("
                        + std::to_string(tstate.getNrow()) + ", " + std::to_string(tstate.getNcol()) + "), expected ("
                        + std::to_string(npoint) + ", " + std::to_string(nstate) + ").";
        throw std::invalid_argument(msg);
    }
    else if (tval.getNcol() != 1) {
        std::string msg = "In computeTransmissionStates(): Invalid size of 'tval'. Received ("
                        + std::to_string(tval.getNrow()) + ", " + std::to_string(tval.getNcol()) + "), expected ("
                        + std::to_string(nstate) + ", 1).";
        throw std::invalid_argument(msg);
    }
    
    // 2. Solve the wave equation and compute the transmission matrix:
    ComplexMatrix tmat(noutputprop, ninputprop);
    transmissionMatrix(tmat);  // We only need the transmission matrix, not the reflection matrix.
    
    // 3. Perform the SVD of the transmission matrix:
    ComplexMatrix u_full(noutputprop, noutputprop), vh_full(ninputprop, ninputprop), vh(nstate, ninputprop);
    RealMatrix sigma(nstate_max, 1);
    svd(tmat, sigma, u_full, vh_full);  // t = U S V^H  -->  t^H t = V S^2 V^H
    
    for (int istate = 0; istate < nstate; istate++) {// Loop on the singular values.
        tval(istate, 0) = sigma(istate, 0)*sigma(istate, 0);  // Compute the transmission eigenvalues from the singular values.
        for (int jmode = 0; jmode < ninputprop; jmode++) {// Loop over the input waveguide modes.
            // Normalize the input state to get approximately unit incident intensity:
            vh(istate, jmode) = vh_full(istate, jmode) * std::conj(I*std::sqrt(2.*ninputprop/(PI*doslattice))*std::sqrt(inputKlh(jmode, 0).real()));
            // Each column of "V" (row of V^H) is the transmission eigenstate (in modal representation).
            // Note that the real part, Re(inputKlh), strips the evanescent modes.
            // Warning: this formula does not take into account multiple input leads.
        }
    }
    
    //std::string filename;
    //uniqueFilename("out/" + sysname + "/vh/vh_", ".ppm", filename);
    //std::cout << TAG_INFO << "Plot output state to file '" << filename << "'...\n";
    //vh.saveImage(filename);
    
    // 4. Compute all the transmission eigenstates:
    computeGreenFunction();  // Ensure the Green function has been computed already (this line does nothing if it is so).
    
    //// Normalize the input state to get approximately unit incident intensity:
    //for (int jmode = 0; jmode < ninputprop; jmode++) {// Loop over the input waveguide modes.
    //    for (int istate = 0; istate < ninputprop; istate++) {// Loop over the input transmission eigenstates.
    //        vh(istate, jmode) *= std::conj(I*std::sqrt(2.*ninputprop/(PI*dosinput))*std::sqrt(inputKlh(jmode, 0).real()));
    //        // Each column of "V" (row of V^H) is the transmission eigenstate (in modal representation).
    //        // Note that the real part, Re(inputKlh), strips the evanescent modes.
    //        // Warning: this formula does not take into account multiple input leads.
    //    }
    //}
    tstate = green * vh.conj();
}

/**
 * Compute several transmission eigenstate profiles associated to given ranges of transmission eigenvalues and add them to the given data matrices.
 * 
 * Arguments:
 * 
 * trange   = Selected ranges of transmission eigenvalues, [Tm-dT, Tm+dT], for which the transmission eigenchannels are desired. Size: (nprofile, 2).
 *            The first column contains the value of the centers of intervals, Tm, and the second column the half-width of the intervals, dT.
 * tprofile = On output, intensity profile (square modulus) of selected transmission eigenstates. Size: (npoint, nprofile).
 *            The transmission eigenstate profiles are normalized so that the incident intensity is equal to 1 (approximately).
 *            Note that, if multiple transmission eigenstates are found in an interval [Tm-dT, Tm+dT], then they are summed up.
 *            If, on the contrary, no eigenstate is found, then the column of "tprofile" is set to zero.
 *            This method adds the profiles directly to "tprofile", so this matrix must be initialized to zero.
 * nsample  = On output, number of found transmission eigenstates in the prescribed intervals. Size: (nprofile, 1).
 *            This method increments the number of found transmission eigenstates, so "nsample" must be initialized to zero.
 * tval     = On output, list of all transmission eigenvalues (ignoring considerations on eigenstates). Size: (min(ninputprop, noutputprop), 1).
 */
void WaveSystem::addITransmission(const RealMatrix& trange, RealMatrix& tprofile, RealMatrix& nsample, RealMatrix& tval) {
    
    // 1. First check for possible errors:
    const int nprofile = trange.getNrow();
    const int ntval = std::min(ninputprop, noutputprop);
    if (trange.getNcol() != 2) {
        throw std::invalid_argument("In addITransmission(): Invalid size of 'trange', expected 2 columns for Tmin and Tmax.");
    }
    else if (tprofile.getNrow() != npoint || tprofile.getNcol() != nprofile) {
        std::string msg = "In addITransmission(): Invalid size of 'tprofile'. Received ("
                        + std::to_string(tprofile.getNrow()) + ", " + std::to_string(tprofile.getNcol()) + "), expected ("
                        + std::to_string(npoint) + ", " + std::to_string(nprofile) + ").";
        throw std::invalid_argument(msg);
    }
    else if (nsample.getNrow() != nprofile || nsample.getNcol() != 1) {
        std::string msg = "In addITransmission(): Invalid size of 'nsample'. Received ("
                        + std::to_string(nsample.getNrow()) + ", " + std::to_string(nsample.getNcol()) + "), expected ("
                        + std::to_string(nprofile) + ", 1).";
        throw std::invalid_argument(msg);
    }
    else if (tval.getNrow() != ntval || tval.getNcol() != 1) {
        std::string msg = "In addITransmission(): Invalid size of 'tval'. Received ("
                        + std::to_string(tval.getNrow()) + ", " + std::to_string(tval.getNcol()) + "), expected ("
                        + std::to_string(ntval) + ", 1).";
        throw std::invalid_argument(msg);
    }
    
    // 2. Compute the singular value decomposition of the transmission matrix:
    ComplexMatrix tmat(noutputprop, ninputprop), u(noutputprop, noutputprop), vh(ninputprop, ninputprop);
    transmissionMatrix(tmat);  // We only need the transmission matrix, not the reflection matrix.
    
    svd(tmat, tval, u, vh);  // Perform the SVD. t = U S V^H  -->  t^H t = V S^2 V^H. 
    // NB: Each row of V^H is the conjugate of a transmission eigenstate (in modal representation). vh(ival, imode)
    
    //double a;
    for (int ival = 0; ival < ntval; ival++) {// Loop on the singular values.
        tval(ival, 0) = tval(ival, 0)*tval(ival, 0);  // Compute the transmission eigenvalues from the singular values.
    }
    
    // 3. Compute the square modulus of selected transmission eigenstates (second most time-consuming operation after the sparse system solution):
    double tm, dt, t;
    dcomplex psi;
    //const dcomplex ampli = 2.*I*std::sqrt(ninputprop/(2*PI*dosinput)); // Old normalization deprecated because depending on the precise value of kh*ninput and thus possibly divergent for resonating input waveguide.
    const dcomplex ampli = 2.*I*std::sqrt(ninputprop/(2*PI*doslattice)); // Normalization of transmission eigenstates such that the incident intensity is 1.
    
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Loop over the desired eigenstate profiles.
        
        tm = trange(iprofile, 0);  // Extract the centroid of the interval [Tm-dT, Tm+dT].
        dt = trange(iprofile, 1);  // Extract the half-width of the interval [Tm-dT, Tm+dT].
        
        for (int ival = 0; ival < ntval; ival++) {// Loop over the transmission eigenvalues to find all eigenstates in the range [Tmin, Tmax].
            
            t = tval(ival, 0);  // Extract the current transmission eigenvalue.
            
            if (std::abs(t - tm) < dt) {// If "t" is in the interval [Tmin, Tmax], then add the profile.
                
                // Compute the square modulus of the transmission state:
                for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
                    
                    // Compute the transmission eigenstate associated to "t":
                    psi = 0.;
                    for (int imode = 0; imode < ninputprop; imode++) {// Loop over the input modes to construct the transmission eigenstate at "ipoint".
                        psi += green(ipoint, imode) * std::conj(vh(ival, imode)) * std::sqrt(inputKlh(imode, 0).real()) * ampli;
                    }
                    tprofile(ipoint, iprofile) += psi.real()*psi.real() + psi.imag()*psi.imag();
                    // Warning: this formula does not take into account multiple input leads.
                }
                nsample(iprofile, 0) += 1;  // Increments the number of found profiles.
            }
        }
    }
}
