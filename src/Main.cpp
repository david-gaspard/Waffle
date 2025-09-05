/****
 * @date Created on 2025-08-13 at 12:43:52 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the main functions of the Waffle program.
 ***/
#include "WaveSystem.hpp"
#include "BaseTools.hpp"
#include <iomanip>

/**
 * @todo TODO LIST:
 * DONE: (1)  In WaveSystem: If possible, compute the exact expression of the free DOS on a square lattice. It is used in setDisorder()...
 * DONE: (2)  In Waffle/Usador: Make the figures for Arthur (see below)...
 * DONE: (3)  In plot_map.py: Fix the partial read bug with the CSV reader which occurs when it is called by the C++ program...
 * DONE: (4)  In WaveSystem: Check in doing math that the normalization of transmission eigenstates is correct.
 * DONE: (5)  Calculate the condition number of the "hamiltonian" for a waveguide by estimating the lowest eigenstate. 
 *            This would indicate whether using iterative methods (instead of direct solvers) is relevant or not...
 * DONE: (6)  In SparseComplexMatrix: Try to optimize the construction of the Hamiltonian.
 *            Check if there is a faster object than "vector" for insertion:
 *            - See "map" (binary tree): faster for insertion than "vector" (no realloc), faster for accessing by key, slower for traversal (no direct linkage).
 *              Probably most appropriate choice if and only if traversal is not too slow compared to "vector".
 *            - See "forward_list" (linked list): faster for insertion than vector (no realloc), faster for traversal, slower for accessing by key 
 *              (no binary search). Both objects should be tested for speed in construction and in conversion for UMFPack.
 * DONE: (7)  Note that if binary trees are indeed faster for sorted insertion/deletion, then "set" would be more appropriate for the points in SquareMesh...
 * DONE: (8)  In Main: Add histogram of transmission eigenvalues. Simply save all raw eigenvalues in a CSV file (each row per disorder realization)...
 * DONE: (9)  Write "plot_histo.py" to plot the histogram of a list of values in given interval [Tmin, Tmax] using Nbins.
 * DONE: (10) In all plot scripts: Extract the header of the CSV file and copy it into the "title" of pgfplots' "axis" environment.
 * DONE: (11) In WaveSystem: Add checkResidual() to verify that the linear system is correctly solved...
 * DONE: (12) In WaveSystem: Maybe remove the evanescent modes from "inputState" and "outputState" because they always give zero transmission eigenvalues...
 * DONE: (13) In SparseComplexMatrix: Implement solveMumps() (no iterative refinement) and compare performance with solveUmfpack()...
 * DONE: (14) Arthur 2025-09-04: Implement the average intensity Ibar(r) without wavefront shaping...
 * DONE: (15) Create script "plot_cut.py" to plot the intensity profile along a cut. A straight line should do the job...
 * DONE: (16) In SparseComplexMatrix: Improve plotImage() to export directly a PNG image instead of a heavy PPM file.
 * 
 * TODO: (17) In SquareMesh: Create addImage(filename) to import PNG in order to facilitate the mesh creation...
 * 
 * (18) Arthur 2025-09-04: For the paper, come back with the double waveguide case (with/without absorber), and compute the eigenstate profile
 *      and the transmission eigenvalue distribution. Test changing the numerical aperture of each of the guides...
 * (19) Arthur 2025-09-04: What happens if we add a non-scattering region in the center of the waveguide ?
 *      We would have D(r) -> infty... What happens in the Usadel equation ?
 *      It would be funny to find a system with transmission eigenstate "opposite" to the average intensity (without shaping). 
 *      Such that maximum of the eigenstate corresponds to the minimum of the average intensity...
 * (20) Do no forget to implement the absorption (holabso != 0) using complex wavenumbers.
 *      Check that the implementation is correct, maybe using the Green function, or the distribution rho(T) (which should match RecurGreen's results)...
 * (21) In Main: Possibly add parallelized taskTransmissionOMP(). Maybe MPI is more relevant, given the memory size of reference 1 Mpx simulations...
 * (22) Arthur 2025-09-04: For the paper, come back on the geometrical interpretation of the Q space and the (theta,eta) parametrization...
 */

/**
 * Defines the overall simulation context.
 */
struct Context {
    WaveSystem sys;
    RealMatrix trange;
};

/***************************************************************************************************
 * SYSTEM CREATION FUNCTIONS
 ***************************************************************************************************/

/**
 * Create a waveguide-shape system of given "length" and "width" (in number of lattice points),
 * and defines the range of transmission values of interest for computing the disorder-averaged profile of transmission eigenchannels. 
 */
Context createWaveguide() {
    
    const int length = 150;   // Number of points in the longitudinal direction, L/h. Better to reach length=width=1000 with kh=1.
    const int width  = 150;   // Number of points in the transverse direction, W/h.
    
    const double dscat = 8.5;  // Scattering depth, L/lscat. Default: dscat=8.5 (in order to get approximately dscat_eff=10).
    const double dabso = 0.;   // Absorption depth, L/labso.
    
    const std::string name = "Waveguide_" + std::to_string(length) + "x" + std::to_string(width) + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, length, 0, width, BND_MIRROR);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, 0, width, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(length, length, 0, width, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/length;
    const double holabso = dabso/length;
    
    RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.980;  trange(0, 1) = 0.020;
    trange(1, 0) = 0.500;  trange(1, 1) = 0.050;
    trange(2, 0) = 0.100;  trange(2, 1) = 0.020;
    trange(3, 0) = 0.005;  trange(3, 1) = 0.001;
    
    return {WaveSystem(name, mesh, kh, holscat, holabso), trange};
}

/**
 * Create a slab system of given length and width.
 */
Context createSlabTransmission1() {
    
    const int length  = 80;   // Number of lattice points in the longitudinal direction, L/h. Default: length=80
    const int width   = 240;  // Number of lattice points in the transverse direction, W/h. Default: width=240
    const int winput  = width/2;  // Width of the input channel. Default: winput=width/2
    const int woutput = width/2;  // Width of the output channel. Default: woutput=width/2
    
    const double dscat = 8.5;  // Scattering depth, L/lscat.
    const double dabso = 0.;  // Absorption depth, L/labso.
    
    const std::string name = "SlabTransmission1_" + std::to_string(length) + "x" + std::to_string(width) 
                           + "_lead" + std::to_string(winput) + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, length, -width/2, width/2, BND_OPEN);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, -winput/2, winput/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(length, length, -woutput/2, woutput/2, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/length;
    const double holabso = dabso/length;
    
    RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.800;  trange(0, 1) = 0.070;
    trange(1, 0) = 0.300;  trange(1, 1) = 0.020;
    trange(2, 0) = 0.100;  trange(2, 1) = 0.010;
    trange(3, 0) = 0.005;  trange(3, 1) = 0.001;
    
    return {WaveSystem(name, mesh, kh, holscat, holabso), trange};
}

/**
 * Create a slab system of given length and width.
 */
Context createSlabTransmission2() {
    
    const int length  = 80;      // Number of lattice points in the longitudinal direction, L/h. Default: length=80
    const int width   = 240;     // Number of lattice points in the transverse direction, W/h. Default: width=240
    const int winput  = width/2; // Width of the input channel. Default: winput=width/2
    const int woutput = width/8; // Width of the output channel. Default: woutput=width/8.
    
    const double dscat = 8.5;  // Scattering depth, L/lscat.
    const double dabso = 0.;  // Absorption depth, L/labso.
    
    const std::string name = "SlabTransmission2_" + std::to_string(length) + "x" + std::to_string(width) 
                           + "_in" + std::to_string(winput) + "_out" + std::to_string(woutput) + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, length, -width/2, width/2, BND_OPEN);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, -winput/2, winput/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(length, length, -woutput/2, woutput/2, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/length;
    const double holabso = dabso/length;
    
    RealMatrix trange(3, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.400;  trange(0, 1) = 0.050;
    trange(1, 0) = 0.100;  trange(1, 1) = 0.010;
    trange(2, 0) = 0.010;  trange(2, 1) = 0.001;
    
    return {WaveSystem(name, mesh, kh, holscat, holabso), trange};
}

/**
 * Create a slab system of given length and width.
 */
Context createSlabTransmission3() {
    
    const int length  = 80;      // Number of lattice points in the longitudinal direction, L/h. Default: length=80
    const int width   = 240;     // Number of lattice points in the transverse direction, W/h. Default: width=240
    const int winput  = width/8; // Width of the input channel. Default: winput=width/2
    const int woutput = width/8; // Width of the output channel. Default: woutput=width/8.
    
    const double dscat = 8.5;  // Scattering depth, L/lscat.
    const double dabso = 0.;  // Absorption depth, L/labso.
    
    const std::string name = "SlabTransmission3_" + std::to_string(length) + "x" + std::to_string(width) 
                           + "_lead" + std::to_string(winput) + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, length, -width/2, width/2, BND_OPEN);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, -winput/2, winput/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(length, length, -woutput/2, woutput/2, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/length;
    const double holabso = dabso/length;
    
    RealMatrix trange(3, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.200;  trange(0, 1) = 0.050;
    trange(1, 0) = 0.100;  trange(1, 1) = 0.010;
    trange(2, 0) = 0.010;  trange(2, 1) = 0.001;
    
    return {WaveSystem(name, mesh, kh, holscat, holabso), trange};
}

/**
 * Create a slab system in remission configuration.
 */
Context createSlabRemission1() {
    
    const int length  = 80;      // Number of lattice points in the longitudinal direction, L/h. Default: length=80
    const int width   = 240;     // Number of lattice points in the transverse direction, W/h. Default: width=240
    const int winput  = width/8; // Width of the input channel. Default: winput=width/2
    const int woutput = width/8; // Width of the output channel. Default: woutput=width/8.
    const int ysep    = width/3; // Separation between the centers of the input and the output.
    
    const double dscat = 8.5;  // Scattering depth, L/lscat.
    const double dabso = 0.;   // Absorption depth, L/labso.
    
    const std::string name = "SlabRemission_" + std::to_string(length) + "x" + std::to_string(width) 
                           + "_lead" + std::to_string(winput) + "_sep" + std::to_string(ysep) + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, length, -width/2, width/2, BND_OPEN);
    
    // Setup boundary conditions:
    mesh.setBoundaryDisk(0, -ysep/2,  winput/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(0, +ysep/2, woutput/2, DIR_WEST, BND_OUTPUT);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/length;
    const double holabso = dabso/length;
    
    RealMatrix trange(3, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.100;  trange(0, 1) = 0.030;
    trange(1, 0) = 0.050;  trange(1, 1) = 0.010;
    trange(2, 0) = 0.010;  trange(2, 1) = 0.002;
    
    return {WaveSystem(name, mesh, kh, holscat, holabso), trange};
}

/**
 * Create a double waveguide with a minute asymmetry.
 */
Context createDoubleWaveguide1() {
    
    const int tlength = 300;          // Number of lattice points in the longitudinal direction, L/h.
    const int twidth  = 300;          // Number of lattice points in the transverse direction, W/h.
    const int xwidth  = tlength/5;    // Width of the corridors near input/output.
    const int ywidth  = twidth/5;     // Width of the corridors of the two waveguides (before the yshift).
    const int yshift  = twidth/50;    // Transverse shift.
    
    const double dscat = 8.5;  // Scattering depth, L/lscat.
    const double dabso = 0.;   // Absorption depth, L/labso.
    
    const std::string name = "DoubleWaveguide1_" + std::to_string(tlength) + "x" + std::to_string(twidth) 
                           + "_cx" + std::to_string(xwidth) + "_cy" + std::to_string(ywidth) + "_yshift" + std::to_string(yshift) 
                           + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, tlength, 0, twidth, BND_MIRROR);
    mesh.removeRectangle(xwidth, tlength-xwidth, ywidth+yshift, twidth-ywidth+yshift);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, 0, twidth, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(tlength, tlength, 0, twidth, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle(xwidth-1, xwidth-1, ywidth+yshift, twidth-ywidth+yshift, DIR_EAST, BND_OPEN);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/tlength;
    const double holabso = dabso/tlength;
    
    RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.880;  trange(0, 1) = 0.010;  //trange(0, 0) = 0.880;  trange(0, 1) = 0.030;
    trange(1, 0) = 0.730;  trange(1, 1) = 0.030;
    trange(2, 0) = 0.500;  trange(2, 1) = 0.020;
    trange(3, 0) = 0.100;  trange(3, 1) = 0.010;
    
    return {WaveSystem(name, mesh, kh, holscat, holabso), trange};
}

/**
 * Create a double waveguide with a minute asymmetry.
 */
Context createDoubleWaveguide2() {
    
    const int tlength = 300;          // Number of lattice points in the longitudinal direction, L/h.
    const int twidth  = 300;          // Number of lattice points in the transverse direction, W/h.
    const int xwidth  = tlength/5;    // Width of the corridors near input/output.
    const int ywidth  = twidth/5;     // Width of the corridors of the two waveguides (before the yshift).
    const int yshift  = twidth/50;    // Transverse shift.
    
    const double dscat = 8.5;  // Scattering depth, L/lscat.
    const double dabso = 0.;   // Absorption depth, L/labso.
    
    const std::string name = "DoubleWaveguide2_" + std::to_string(tlength) + "x" + std::to_string(twidth) 
                           + "_cx" + std::to_string(xwidth) + "_cy" + std::to_string(ywidth) + "_yshift" + std::to_string(yshift) 
                           + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, tlength, 0, twidth, BND_MIRROR);
    mesh.removeRectangle(xwidth, tlength-xwidth, ywidth+yshift, twidth-ywidth+yshift);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, 0, twidth, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(tlength, tlength, 0, twidth, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/tlength;
    const double holabso = dabso/tlength;
    
    RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.980;  trange(0, 1) = 0.020;
    trange(1, 0) = 0.880;  trange(1, 1) = 0.030;
    trange(2, 0) = 0.730;  trange(2, 1) = 0.030;
    trange(3, 0) = 0.500;  trange(3, 1) = 0.020;
    
    return {WaveSystem(name, mesh, kh, holscat, holabso), trange};
}

/***************************************************************************************************
 * COMPUTATION OF TRANSMISSION EIGENSTATES
 ***************************************************************************************************/

/**
 * Normalize the disorder averages of transmission eigenstates profile contained in "tprofile" knowing the corresponding number of samples "nsample".
 */
void normalizeTProfile(RealMatrix& tprofile, const RealMatrix& nsample) {
    
    // 1. Check for possible errors:
    const int npoint = tprofile.getNrow();
    const int nprofile = tprofile.getNcol();
    if (nsample.getNrow() != nprofile || nsample.getNcol() != 1) {
        std::string msg = "In normalizeTProfile(): Invalid dimensions of nsample, received (" 
                        + std::to_string(nsample.getNrow()) + ", " + std::to_string(nsample.getNcol()) + "), expected (" 
                        + std::to_string(nprofile) + ", 1).";
        throw std::invalid_argument(msg);
    }
    
    // 2. Compute the average profiles by dividing by the number of samples:
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Loop over the transmission eigenstate profiles.
        if (nsample(iprofile, 0) > MEPS) {// Only normalize the profile if the number of samples is nonzero.
            for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
                tprofile(ipoint, iprofile) /= nsample(iprofile, 0);
            }
        }
    }
}

/**
 * Save the transmission profile to a CSV file and call external scripts to plot it.
 */
void plotTProfile(const RealMatrix& tprofile, const RealMatrix& trange, const RealMatrix& nsample,
                  const RealMatrix& tvalstore, const WaveSystem& sys, const double ctime) {
    
    // 1. First compute some metadata:
    const int nseed = tvalstore.getNrow();
    const int nprofile = tprofile.getNcol();
    const int ntval = std::min(sys.getNInputProp(), sys.getNOutputProp());  // Get the number of propagationg modes (nonzero transmission eigenvalues).
    const double tavg = tvalstore.sum()/(ntval*static_cast<double>(nseed));  // Compute the average transmission probability.
    const double dscat_eff = (PI/2.) * (1./tavg - 1.);   // Deduce the effective scattering depth (as if in a waveguide).
    
    // 2. Save the transmission eigenstate profiles to a file and plot them:
    std::string info = "Transmission eigenstate profiles for Trange=[";
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Save the ranges of transmission eigenvalues for which 
        info += " " + std::to_string(trange(iprofile, 0)) + "+-" + std::to_string(trange(iprofile, 1)) + " ";
    }
    info += "], found Nsample=[";
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Save the ranges of transmission eigenvalues for which 
        info += " " + std::to_string(nsample(iprofile, 0)) + " ";
    }
    info += "], Nseed=" + std::to_string(nseed) + ", Tavg=" + std::to_string(tavg) + ", L/lscat_eff=" + std::to_string(dscat_eff) + ", Computation_time=" + std::to_string(ctime) + " s.";
    
    sys.plotIntensity(tprofile, info, "tprofile");  // Plot the average transmission eigenstate profiles.
}

/**
 * Save the given samples to a file and call external scripts to plot the normalized histogram.
 * each row in "tvalstore" are the eigenvalue samples corresponding to one realization of the disorder
 */
void plotHistogram(const RealMatrix& tvalstore, const WaveSystem& sys, const double ctime) {
    
    const int nseed = tvalstore.getNrow();  // Number of realizations of the disorder.
    const int ntval = tvalstore.getNcol();  // Number of nonzero eigenvalues (also the number of propagating modes).
    const double tavg = tvalstore.sum()/(ntval*static_cast<double>(nseed));  // Compute the average transmission probability.
    const double dscat_eff = (PI/2.) * (1./tavg - 1.);   // Deduce the effective scattering depth (as if in a waveguide).
    const char* sep = ", ";  // Separator used between entries of the CSV file.
    const int prec = 16;     // Precision used in printing double precision values.
    const std::string filename = sys.uniqueFile("tspectrum", ".csv");
    const double fsize = (prec+4.) * ntval * static_cast<double>(nseed);  // Roughly estimated file size in bytes (octets).
    std::cout << TAG_INFO << "Save samples to file '" << filename << "', size ~" << (fsize/1e6) << " Mo...\n";
    
    std::ofstream ofs(filename);    // Open the output file.
    ofs << std::setprecision(prec); // Set the printing precision.
    writeTimestamp(ofs, "%% ");     // Apply a timestamp at the beginning.
    
    for (const std::string& line : sys.summary()) {// Write the summary to the file header.
        ofs << "%% " << line << "\n";
    }
    std::string info = "Transmission eigenvalues for Nprop=" + std::to_string(ntval) + ", Nseed=" + std::to_string(nseed) 
        + ", Tavg=" + std::to_string(tavg) + ", L/lscat_eff=" + std::to_string(dscat_eff) + ", Computation_time=" + std::to_string(ctime) + " s.";
    
    ofs << "%% Info: " << info << "\nT_0";
    
    for (int ival = 1; ival < ntval; ival++) {// Loop over the input modes to finish the column names.
        ofs << sep << "T_" << ival;
    }
    ofs << "\n";
    
    for (int iseed = 0; iseed < nseed; iseed++) {// Loop over the rows of tvalstore.
        ofs << tvalstore(iseed, 0);
        for (int ival = 1; ival < ntval; ival++) {// Loop over the samples.
            ofs << sep << tvalstore(iseed, ival);
        }
        ofs << "\n";
    }
    ofs.close();  // Close the stream before calling an external script (this may cause I/O trouble).
    
    std::string cmd("plot/plot_histo.py " + filename);
    std::cout << TAG_EXEC << cmd << "\n";
    if (std::system(cmd.c_str())) {
        std::cout << TAG_WARN << "The plot script returned an error.\n";
    }
}

/**
 * Compute all quantities related to transmission for the given WaveSystem "sys", in particular the transmission-eigenvalue distribution,
 * and the intensity profile of selected transmission eigenchannels (aka eigenstates).
 * Save the results to files with automated names, and call external plot scripts.
 */
void taskTransmissionSerial(WaveSystem& sys, RealMatrix& trange, const int nseed) {
    
    const int npoint = sys.getNPoint();
    const int nprofile = trange.getNrow();  // Get the desired number of profiles.
    const int ntval = std::min(sys.getNInputProp(), sys.getNOutputProp());  // Get the number of nonzero transmission eigenvalues (propagating modes).
    RealMatrix tprofile(npoint, nprofile), nsample(nprofile, 1), tval(ntval, 1), tvalstore(nseed, ntval);
    
    const std::string msg = "Tprofile, serial";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
        
        sys.setDisorder(iseed+1); // Avoid seed=0 for safety.
        sys.addTransmissionProfiles(trange, tprofile, nsample, tval);
        
        for (int ival = 0; ival < ntval; ival++) {// Save the nonzero transmission eigenvalues in the matrix tvalstore.
            tvalstore(iseed, ival) = tval(ival, 0); // Copy the transmission eigenvalues.
        }
        
        printProgressBar(iseed+1, nseed, msg, start);
        
        //if (nsample(0, 0) >= 1) break;  // Conditional break (useful to find exactly one eigenstate).
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save data to files and plot:
    normalizeTProfile(tprofile, nsample);  // Normalize the averages of transmission eigenstate profiles.
    plotTProfile(tprofile, trange, nsample, tvalstore, sys, ctime);
    plotHistogram(tvalstore, sys, ctime);
}

/***************************************************************************************************
 * COMPUTATION OF MODE-AVERAGED INTENSITY
 ***************************************************************************************************/

/**
 * Compute the mode-average intensity averaged over "nseed" realizations of the disorder.
 * Save the results to files with automated names, and call external plot scripts.
 */
void taskAverageIntensitySerial(WaveSystem& sys, const int nseed) {
    
    const int npoint = sys.getNPoint();
    RealMatrix ibar(npoint, 1);
    
    const std::string msg = "Ibar, serial";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
        
        sys.setDisorder(iseed+1); // Avoid seed=0 for safety.
        sys.addAverageIntensity(ibar);
        printProgressBar(iseed+1, nseed, msg, start);
    }
    
    // Normalize the disorder average:
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
        ibar(ipoint, 0) /= nseed;
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save the data and plot:
    const std::string info = "Mode-averaged intensity Ibar for Nseed=" + std::to_string(nseed) + ", Computation_time=" + std::to_string(ctime) + " s.";
    sys.plotIntensity(ibar, info, "ibar");
}

/***************************************************************************************************
 * MAIN FUNCTION OF THE PROGRAM
 ***************************************************************************************************/
int main(int argc, char** argv) {
    
    std::cout << "****** This is " << PROGRAM_COPYRIGHT << " ******\n";
    
    // Create the simulation context:
    Context ctx = createWaveguide();
    //Context ctx = createSlabTransmission1();
    //Context ctx = createSlabTransmission2();
    //Context ctx = createSlabTransmission3();
    //Context ctx = createSlabRemission1();
    //Context ctx = createDoubleWaveguide1();
    //Context ctx = createDoubleWaveguide2();
    
    ctx.sys.printSummary();
    ctx.sys.infoHamiltonian();
    
    // Checking operations:
    //ctx.sys.setDisorder(1);
    //ctx.sys.checkResidual();
    //ctx.sys.checkUnitarity(true);
    //ctx.sys.plotMesh();
    //ctx.sys.plotInputState();
    //ctx.sys.plotOutputState();
    //ctx.sys.plotHamiltonian();
    //ctx.sys.plotGreenFunction();
    //ctx.sys.plotTransmissionStates();
    
    //const int nseed = 100; // Number of random realizations of the disorder used for averaging. Recommended for high quality: 10^4.
    
    //taskTransmissionSerial(ctx.sys, ctx.trange, nseed);
    //taskTransmissionOMP(ctx.sys, ctx.trange, nseed);
    //taskAverageIntensitySerial(ctx.sys, nseed);
    
    return 0;
}
