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
 * DONE: (17) In SquareMesh: Create addImage(filename) to import PNG in order to facilitate the mesh creation...
 * DONE: (18) Arthur 2025-09-04: For the paper, come back with the double waveguide case (with/without absorber), and compute the eigenstate profile
 *       and the transmission eigenvalue distribution. Test changing the numerical aperture of each of the guides...
 * DONE: (19) In Main: Possibly add parallelized taskTransmissionOMP(). Maybe MPI is more relevant, given the memory size of reference 1 Mpx simulations...
 * DONE: (20) Modify plot_histo.py, plot_cut.py, plot_proj.py to output CSV files instead of writing data directly into TikZ files.
 *            This will allow the simulations parameters in headers to be preserved as locally as possible.
 * DONE: (21) Arthur 2025-09-11: Implement average intensity for a plane wave input in order to show the comparison in the paper.
 *            In the code, make clear the distinction between Lambertian input and plane wave input.
 * DONE: (22) In Usador: Add computation time, especially for the computation of the Q field at specific T.
 * ----------------------------------------------------------------------
 * TODO (2025-09-25): 
 * The normalization problem of 2025-09-24 is now considered as solved.
 * However, some safety measures and revision should be done to avoid encountering this problem again:
 * (1) Create a safety system to avoid that the input waveguide resonates :
 *      kh_fix = 2*sin( pi/(2*(Ny+1)) * ( floor( (2*(Ny + 1))/pi * arcsin(kh/2) ) + 0.3 ) );
 * This should give approximately: dosinput = dosfree.
 * (2) Remove estimate of dscat_exact (with accurate extrapolation length because Usador does not use this value (this can create confusion).
 * (3) Improve estimate of dscat_approx (with approximation extrapolation length) with the computation of statistical uncertainties (as in RecurGreen).
 * ----------------------------------------------------------------------
 * (23) Perform benchmark simulations (+ comparison with Usadel) for :
 *      (A) Square waveguide 500x500, 
 *      (B) Slab transmission 300x900, equal size input and output,
 *      (C) Slab transmission 300x900, large input, small output,
 *      (D) Slab remission 300x900, 
 *      (E) Another geometry, maybe the double waveguide, the maze, the "random fiber", or simply the random cavity (to be discussed)...
 * 
 * (24) Arthur 2025-09-17: Consider the trinity of waveguides with absorbers: Tmax, Rmin, Abso_max (=Smin).
 *      To this end, we need (1) to create the "inout" boundary condition, (2) to implement local absorption. And do the same in Usador.
 * (25) Do no forget to implement the absorption (holabso != 0) using complex wavenumbers.
 *      Check that the implementation is correct, maybe using the Green function, or the distribution rho(T) (which should match RecurGreen's results)...
 * (26) Arthur 2025-09-04: What happens if we add a non-scattering region in the center of the waveguide ?
 *      We would have D(r) -> infty... What happens in the Usadel equation ?
 *      It would be funny to find a system with transmission eigenstate "opposite" to the average intensity (without shaping). 
 *      Such that maximum of the eigenstate corresponds to the minimum of the average intensity...
 * (27) Arthur 2025-09-04: For the paper, come back on the geometrical interpretation of the Q space and the (theta,eta) parametrization...
 */

/**
 * Defines the overall simulation context.
 */
struct Context {
    WaveSystem sys;
    RealMatrix trange;
};

/***************************************************************************************************
 * TEST ROUTINES
 ***************************************************************************************************/

/**
 * Check the unitarity of the propagation for different realizations of the disorder in parallel.
 * Thisfunction mainly serves as a parallelization test, especially to check the compatibility with MUMPS.
 */
void taskCheckUnitarityOMP(const WaveSystem& sys, const int nseed, const int seed0, const int nthread) {
    
    const bool showtval = false; // Show the transmission eigenvalues or not.
    const std::string msg = "checkUnitarity, " + std::to_string(nthread) + " thr";
    int cjob = 0;
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    #pragma omp parallel num_threads(nthread)
    {
        WaveSystem sys_loc(sys); // Deep copy of the system on each thread (OMP private).
        
        #pragma omp for schedule(static)
        for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
            
            sys_loc.setDisorder(seed0 + iseed); // Avoid seed=0 for safety.
            sys_loc.checkUnitarity(showtval);
            
            // Critical section to print the progress bar:
            #pragma omp critical
            {
                cjob++;
                printProgressBar(cjob, nseed, msg, start);
            }
        }
    }
    endProgressBar(start);
}

/***************************************************************************************************
 * COMPUTATION OF TRANSMISSION EIGENSTATES
 ***************************************************************************************************/

/**
 * Normalize the disorder averages of transmission eigenstates profile contained in "tprofile" knowing the corresponding number of samples "nsample".
 */
void normalizeITransmission(RealMatrix& tprofile, const RealMatrix& nsample) {
    
    // 1. Check for possible errors:
    const int npoint = tprofile.getNrow();
    const int nprofile = tprofile.getNcol();
    if (nsample.getNrow() != nprofile || nsample.getNcol() != 1) {
        std::string msg = "In normalizeITransmission(): Invalid dimensions of nsample, received (" 
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
 * Compute the average transmission probability (with uncertainty) and estimate the effective scattering thickness (with uncertainty)
 * under the assumption that the system has a waveguide geometry.
 */
void estimateScatteringDepth(const RealMatrix& tvalstore, double& tavg, double& tuct, double& dscat_eff, double& dscat_uct) {
    
    const int nseed = tvalstore.getNrow();
    const int ntval = tvalstore.getNcol();
    RealMatrix tavgstore(nseed, 1);
    
    for (int iseed = 0; iseed < nseed; iseed++) {// Loop over the realizations of the disorder.
        tavgstore(iseed, 0) = 0.;
        for (int ival = 0; ival < ntval; ival++) {
            tavgstore(iseed, 0) += tvalstore(iseed, ival);
        }
        tavgstore(iseed, 0) /= ntval;
    }
    
    tavg = tavgstore.mean();    // Compute the average transmission probability.
    tuct = tavgstore.stddev();  // Uncertainty over the average value "tavg".
    dscat_eff = 2*EXTRAPOLEN_APPROX * (1./tavg - 1.);   // Deduce the effective scattering depth (as if in a waveguide).
    dscat_uct = 2*EXTRAPOLEN_APPROX * tuct/(tavg*tavg);  // Deduce the uncertainty over the scattering depth.
    //const double dscat_eff_exact  = 2*EXTRAPOLEN_EXACT  * (1./tavg - 1.);   // Deduce the effective scattering depth using ht eexact extrapolation length.
    // Using the exact version is not recommended because we need to use the same formula in every program to be consistent and take advantage of error cancellations on the extrpolation length. This is especially when comparing with the Usadel equation.
}

/**
 * Save the intensity profile of transmission eigenstate to a CSV file and call external scripts to plot it.
 */
void plotITransmission(const RealMatrix& tprofile, const RealMatrix& trange, const RealMatrix& nsample, const RealMatrix& tvalstore, 
                       const WaveSystem& sys, const int nseed, const int seed0, const int nthread, const double ctime) {
    
    if (tvalstore.getNrow() != nseed) {// Check for possible errors.
        std::string msg = "In plotITransmission(): Invalid number of seeds, received tvalstore.nrow=" 
                        + std::to_string(tvalstore.getNrow()) + ", but nseed=" + std::to_string(nseed);
        throw std::invalid_argument(msg);
    }
    
    const int nprofile = tprofile.getNcol();
    double tavg, tuct, dscat_eff, dscat_uct;
    estimateScatteringDepth(tvalstore, tavg, tuct, dscat_eff, dscat_uct);
    
    std::string info = "Transmission eigenstate profiles for Trange=[";
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Save the ranges of transmission eigenvalues for which 
        info += " " + std::to_string(trange(iprofile, 0)) + "+-" + std::to_string(trange(iprofile, 1)) + " ";
    }
    info += "]\n%% Found Nsample=[";
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Save the ranges of transmission eigenvalues for which 
        info += " " + std::to_string(static_cast<int>(nsample(iprofile, 0))) + " ";
    }
    info += "], Nseed=" + std::to_string(nseed) + ", Seed0=" + std::to_string(seed0) + ", Tavg=" + std::to_string(tavg) + "+-" + std::to_string(tuct)
        + ", L/lscat_eff=" + std::to_string(dscat_eff) + "+-" + std::to_string(dscat_uct)
        + ", Computation_time=" + std::to_string(ctime) + " s, Nthread=" + std::to_string(nthread) + ".";
    
    sys.plotIntensity(tprofile, info, "itprofile");  // Plot the average transmission eigenstate profiles.
}

/**
 * Save the given samples to a file and call external scripts to plot the normalized histogram.
 * each row in "tvalstore" are the eigenvalue samples corresponding to one realization of the disorder
 */
void plotHistogram(const RealMatrix& tvalstore, const WaveSystem& sys, const int nseed, const int seed0, const int nthread, const double ctime) {
    
    if (tvalstore.getNrow() != nseed) {// Check for possible errors.
        std::string msg = "In plotITransmission(): Invalid number of seeds, received tvalstore.nrow=" 
                        + std::to_string(tvalstore.getNrow()) + ", but nseed=" + std::to_string(nseed);
        throw std::invalid_argument(msg);
    }
    
    const int ntval = tvalstore.getNcol();  // Number of nonzero eigenvalues (also the number of propagating modes).
    double tavg, tuct, dscat_eff, dscat_uct;
    estimateScatteringDepth(tvalstore, tavg, tuct, dscat_eff, dscat_uct);
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
        + ", Seed0=" + std::to_string(seed0) + ", Tavg=" + std::to_string(tavg) + "+-" + std::to_string(tuct)
        + ", L/lscat_eff=" + std::to_string(dscat_eff) + "+-" + std::to_string(dscat_uct)
        + ", Computation_time=" + std::to_string(ctime) + " s, Nthread=" + std::to_string(nthread) + ".";
    
    ofs << "%% Info: " << info << "\nT0";
    
    for (int ival = 1; ival < ntval; ival++) {// Loop over the input modes to finish the column names.
        ofs << sep << "T" << ival;
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
void taskITransmissionSerial(WaveSystem& sys, RealMatrix& trange, const int nseed, const int seed0) {
    
    const int npoint = sys.getNPoint();
    const int nprofile = trange.getNrow();  // Get the desired number of profiles.
    const int ntval = std::min(sys.getNInputProp(), sys.getNOutputProp());  // Get the number of nonzero transmission eigenvalues (propagating modes).
    RealMatrix tprofile(npoint, nprofile), nsample(nprofile, 1), tval(ntval, 1), tvalstore(nseed, ntval);
    
    const std::string msg = "Iavg_transmission, serial";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
        
        sys.setDisorder(seed0 + iseed); // Avoid seed=0 for safety.
        sys.addITransmission(trange, tprofile, nsample, tval);
        
        for (int ival = 0; ival < ntval; ival++) {// Save the nonzero transmission eigenvalues in the matrix tvalstore.
            tvalstore(iseed, ival) = tval(ival, 0); // Copy the transmission eigenvalues.
        }
        
        printProgressBar(iseed+1, nseed, msg, start);
        
        //if (nsample(0, 0) >= 1) break;  // Conditional break (useful to find exactly one eigenstate).
    }
    
    // Normalize the averages of transmission eigenstate profiles:
    normalizeITransmission(tprofile, nsample);
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save data to files and plot:
    plotITransmission(tprofile, trange, nsample, tvalstore, sys, nseed, seed0, 1, ctime);
    plotHistogram(tvalstore, sys, nseed, seed0, 1, ctime);
}

/**
 * Compute all quantities related to transmission for the given WaveSystem "sys", in particular the transmission-eigenvalue distribution,
 * and the intensity profile of selected transmission eigenchannels (aka eigenstates).
 * Save the results to files with automated names, and call external plot scripts.
 * This version uses multithreading with OpenMP to perform the averaging over the disorder.
 */
void taskITransmissionOMP(WaveSystem& sys, RealMatrix& trange, const int nseed, const int seed0, const int nthread) {
    
    const int npoint = sys.getNPoint();
    const int nprofile = trange.getNrow();  // Get the desired number of profiles.
    const int ntval = std::min(sys.getNInputProp(), sys.getNOutputProp());  // Get the number of nonzero transmission eigenvalues (propagating modes).
    RealMatrix tprofile(npoint, nprofile), nsample(nprofile, 1), tval(ntval, 1), tvalstore(nseed, ntval);
    
    int cjob = 0; // Current number of completed jobs (i.e., realizations of the disorder).
    const std::string msg = "Iavg_transmission, " + std::to_string(nthread) + " thr";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    #pragma omp parallel num_threads(nthread)
    {
        WaveSystem sys_loc(sys); // Deep copy of the system on each thread (OMP private).
        RealMatrix tprofile_loc(tprofile), nsample_loc(nsample), tval_loc(tval); // Local data (OMP private).
        
        #pragma omp for schedule(dynamic)
        for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
            
            sys_loc.setDisorder(seed0 + iseed); // Avoid seed=0 for safety.
            sys_loc.addITransmission(trange, tprofile_loc, nsample_loc, tval_loc);
            
            for (int ival = 0; ival < ntval; ival++) {// Save the nonzero transmission eigenvalues in the matrix tvalstore.
                tvalstore(iseed, ival) = tval_loc(ival, 0); // Copy the transmission eigenvalues.
            }
            
            // Critical section to print the progress bar:
            #pragma omp critical
            {
                cjob++;
                printProgressBar(cjob, nseed, msg, start);
            }
        }
        
        // Critical section to gather all data together:
        #pragma omp critical
        {
            tprofile += tprofile_loc;
            nsample += nsample_loc;
        }
    }
    
    // Normalize the averages of transmission eigenstate profiles:
    normalizeITransmission(tprofile, nsample);
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save data to files and plot:
    plotITransmission(tprofile, trange, nsample, tvalstore, sys, nseed, seed0, nthread, ctime);
    plotHistogram(tvalstore, sys, nseed, seed0, nthread, ctime);
}

/***************************************************************************************************
 * COMPUTATION OF AVERAGED INTENSITY
 ***************************************************************************************************/

/**
 * Compute the intensity corresponding to an isotropic (Lambertian) input state averaged over "nseed" realizations of the disorder.
 * Save the results to files with automated names, and call external plot scripts.
 */
void taskIIsotropicSerial(WaveSystem& sys, const int nseed, const int seed0) {
    
    const int npoint = sys.getNPoint();
    RealMatrix iavgiso(npoint, 1);
    
    const std::string msg = "Iavg_isotropic, serial";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
        
        sys.setDisorder(seed0 + iseed); // Avoid seed=0 for safety.
        sys.addIIsotropic(iavgiso);
        printProgressBar(iseed+1, nseed, msg, start);
    }
    
    // Normalize the disorder average:
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
        iavgiso(ipoint, 0) /= nseed;
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save the data and plot:
    const std::string info = "Intensity for isotropic input with Nseed=" + std::to_string(nseed) 
        + ", Seed0=" + std::to_string(seed0) + ", Computation_time=" + std::to_string(ctime) + " s (serial).";
    sys.plotIntensity(iavgiso, info, "iavg_iso");
}

/**
 * Compute the intensity corresponding to an isotropic (Lambertian) input state averaged over "nseed" realizations of the disorder.
 * Save the results to files with automated names, and call external plot scripts.
 * This version uses multithreading with OpenMP.
 */
void taskIIsotropicOMP(WaveSystem& sys, const int nseed, const int seed0, const int nthread) {
    
    const int npoint = sys.getNPoint();
    RealMatrix iavgiso(npoint, 1);
    
    int cjob = 0; // Current number of completed jobs (i.e., realizations of the disorder).
    const std::string msg = "Iavg_isotropic, " + std::to_string(nthread) + " thr";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    #pragma omp parallel num_threads(nthread)
    {
        WaveSystem sys_loc(sys); // Deep copy of the system on each thread (OMP private).
        RealMatrix iavgiso_loc(npoint, 1);
        
        #pragma omp for schedule(dynamic)
        for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
            
            sys_loc.setDisorder(seed0 + iseed); // Avoid seed=0 for safety.
            sys_loc.addIIsotropic(iavgiso_loc);
            
            // Critical section to print the progress bar:
            #pragma omp critical
            {
                cjob++;
                printProgressBar(cjob, nseed, msg, start);
            }
        }
        
        // Critical section to gather all data together:
        #pragma omp critical
        {
            iavgiso += iavgiso_loc;
        }
    }
    
    // Normalize the disorder average:
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
        iavgiso(ipoint, 0) /= nseed;
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save the data and plot:
    const std::string info = "Intensity for isotropic input with Nseed=" + std::to_string(nseed) 
        + ", Seed0=" + std::to_string(seed0) + ", Computation_time=" + std::to_string(ctime) + " s, Nthread=" + std::to_string(nthread) + ".";
    sys.plotIntensity(iavgiso, info, "iavg_iso");
}

/**
 * Compute the intensity corresponding to an incident input mode averaged over "nseed" realizations of the disorder.
 * Save the results to files with automated names, and call external plot scripts.
 * This version uses multithreading with OpenMP.
 */
void taskIModeOMP(WaveSystem& sys, const int imode, const int nseed, const int seed0, const int nthread) {
    
    const int npoint = sys.getNPoint();
    RealMatrix iavgmode(npoint, 1);
    
    int cjob = 0; // Current number of completed jobs (i.e., realizations of the disorder).
    const std::string msg = "Iavg_mode, " + std::to_string(nthread) + " thr";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    #pragma omp parallel num_threads(nthread)
    {
        WaveSystem sys_loc(sys); // Deep copy of the system on each thread (OMP private).
        RealMatrix iavgmode_loc(npoint, 1);
        
        #pragma omp for schedule(dynamic)
        for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
            
            sys_loc.setDisorder(seed0 + iseed); // Avoid seed=0 for safety.
            sys_loc.addIMode(iavgmode_loc, imode);
            
            // Critical section to print the progress bar:
            #pragma omp critical
            {
                cjob++;
                printProgressBar(cjob, nseed, msg, start);
            }
        }
        
        // Critical section to gather all data together:
        #pragma omp critical
        {
            iavgmode += iavgmode_loc;
        }
    }
    
    // Normalize the disorder average:
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
        iavgmode(ipoint, 0) /= nseed;
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save the data and plot:
    const std::string info = "Intensity for plane input with Imode=" + std::to_string(imode) + ", Nseed=" + std::to_string(nseed) 
        + ", Seed0=" + std::to_string(seed0) + ", Computation_time=" + std::to_string(ctime) + " s, Nthread=" + std::to_string(nthread) + ".";
    sys.plotIntensity(iavgmode, info, "iavg_mode");
}

/***************************************************************************************************
 * COMPUTATION OF ALL INSTANCES
 ***************************************************************************************************/

/**
 * Compute all the field intensities and other instances in a single run using multithreading with OpenMP.
 */
void taskIAllOMP(WaveSystem& sys, const RealMatrix& trange, const int imode, const int nseed, const int seed0, const int nthread) {
    
    const int npoint = sys.getNPoint();
    const int nprofile = trange.getNrow();  // Get the desired number of profiles.
    const int ntval = std::min(sys.getNInputProp(), sys.getNOutputProp());  // Get the number of nonzero transmission eigenvalues (propagating modes).
    RealMatrix tprofile(npoint, nprofile), nsample(nprofile, 1), tval(ntval, 1), tvalstore(nseed, ntval);
    RealMatrix iavgiso(npoint, 1), iavgmode(npoint, 1);
    
    int cjob = 0; // Current number of completed jobs (i.e., realizations of the disorder).
    std::string info = "Iall, " + std::to_string(nthread) + " thr";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    #pragma omp parallel num_threads(nthread)
    {
        WaveSystem sys_loc(sys); // Deep copy of the system on each thread (OMP private).
        RealMatrix tprofile_loc(tprofile), nsample_loc(nsample), tval_loc(tval); // Local data (OMP private).
        RealMatrix iavgiso_loc(npoint, 1), iavgmode_loc(npoint, 1); // Local data (OMP private).
        
        #pragma omp for schedule(dynamic)
        for (int iseed = 0; iseed < nseed; iseed++) {// Loop over realizations of the disorder.
            
            sys_loc.setDisorder(seed0 + iseed);
            sys_loc.addITransmission(trange, tprofile_loc, nsample_loc, tval_loc);
            sys_loc.addIIsotropic(iavgiso_loc);
            sys_loc.addIMode(iavgmode_loc, imode);
            
            for (int ival = 0; ival < ntval; ival++) {// Save the nonzero transmission eigenvalues in the matrix tvalstore.
                tvalstore(iseed, ival) = tval_loc(ival, 0); // Copy the transmission eigenvalues.
            }
            
            // Critical section to print the progress bar:
            #pragma omp critical
            {
                cjob++;
                printProgressBar(cjob, nseed, info, start);
            }
        }
        
        // Critical section to gather all data together:
        #pragma omp critical
        {
            tprofile += tprofile_loc;
            nsample += nsample_loc;
            iavgiso += iavgiso_loc;
            iavgmode += iavgmode_loc;
        }
    }
    
    // Normalize the disorder average:
    normalizeITransmission(tprofile, nsample);  // Normalize the averages of transmission eigenstate profiles.
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
        iavgiso(ipoint, 0) /= nseed;
        iavgmode(ipoint, 0) /= nseed;
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    // Save data to files and plot:
    plotITransmission(tprofile, trange, nsample, tvalstore, sys, nseed, seed0, nthread, ctime);
    plotHistogram(tvalstore, sys, nseed, seed0, nthread, ctime);
    info = "Intensity for isotropic input with Nseed=" + std::to_string(nseed) + ", Seed0=" + std::to_string(seed0) 
         + ", Computation_time=" + std::to_string(ctime) + " s, Nthread=" + std::to_string(nthread) + ".";
    sys.plotIntensity(iavgiso, info, "iavg_iso");
    info = "Intensity for plane input with Imode=" + std::to_string(imode) + ", Nseed=" + std::to_string(nseed) + ", Seed0=" + std::to_string(seed0) 
         + ", Computation_time=" + std::to_string(ctime) + " s, Nthread=" + std::to_string(nthread) + ".";
    sys.plotIntensity(iavgmode, info, "iavg_mode");
}

/***************************************************************************************************
 * MAIN FUNCTION OF THE PROGRAM
 ***************************************************************************************************/
int main(int argc, char** argv) {
    
    std::cout << "****** This is " << PROGRAM_COPYRIGHT << " ******\n";
    
    // Create the system from a PNG image:
    //SquareMesh mesh("model/waveguide_30x20.png");
    //SquareMesh mesh("model/waveguide_302x300.png"); // Currently standard waveguide.
    //SquareMesh mesh("model/waveguide_602x600.png");
    //SquareMesh mesh("model/waveguide_1202x1200.png");
    //SquareMesh mesh("model/double-guide-abso-sym_642x384.png");
    //SquareMesh mesh("model/double-guide-abso-shift-15_642x384.png");
    //SquareMesh mesh("model/maze_706x513.png");
    //SquareMesh mesh("model/slab-transmission-1_299x893.png");
    //SquareMesh mesh("model/ring-guide_700x700.png");
    //SquareMesh mesh("model/arched-guide_1002x750.png");
    //SquareMesh mesh("model/arched-guide-absorber_1002x750.png");
    //SquareMesh mesh("model/slab-guide_152x450.png"); // Used for L/lscat calibration.
    //SquareMesh mesh("model/slab-guide_302x900.png"); // Used for L/lscat calibration.
    //SquareMesh mesh("model/slab-guide_602x600.png"); // Used for L/lscat calibration.
    //SquareMesh mesh("model/slab-guide_101x297.png");
    //SquareMesh mesh("model/slab-transmission-1_302x902.png");
    //SquareMesh mesh("model/slab-transmission-2_302x902.png");
    //SquareMesh mesh("model/slab-transmission-focus_302x902.png");
    //SquareMesh mesh("model/slab-tm-ar3-in7-out1_302x902.png");
    //SquareMesh mesh("model/slab-tm-ar3-in7-out7_302x902.png");
    //SquareMesh mesh("model/slab-tm-ar3-in7-out15_302x902.png");
    //SquareMesh mesh("model/slab-remission_302x902.png");
    //SquareMesh mesh("model/slab-guide_102x300.png");
    //SquareMesh mesh("model/slab-tm-ar3-in1-out1_102x302.png"); // Test the exponential decay of the profile in the transverse direction.
    //SquareMesh mesh("model/slab-tm-ar9-in1-out1_102x902.png"); // Test the exponential decay of the profile in the transverse direction.
    //SquareMesh mesh("model/waveguide-obstacle_302x300.png");
    //SquareMesh mesh("model/waveguide_502x500.png");
    //SquareMesh mesh("model/slab-remission_152x452.png");
    //SquareMesh mesh("model/slab-tm-ar3-div15-in45-out5_152x452.png");
    //SquareMesh mesh("model/slab-tm-ar3-div15-in5-out5_302x902.png");  // Currently standard slab (for dscat_eff=20, dscat=16.40).
    //SquareMesh mesh("model/slab-tm-ar3-div15-in9-out9_302x902.png");
    //SquareMesh mesh("model/shrinking-guide_602x1000.png");
    //SquareMesh mesh("model/random-cavity-4-leads-1_600x385.png");
    //SquareMesh mesh("model/waveguide-obstacle-sym_302x300.png");
    //SquareMesh mesh("model/year-2026_1762x578.png");
    SquareMesh mesh("model/constriction-1c4_1000x1000.png");
    
    /**
     * TABLE of effective scattering thickness for 300x900 slab waveguides with kh=1 and holscat=dscat/300:
     * dscat    dscat_eff_observed | OLD FIT = 1.13763 dscat + 0.00122407 dscat^2 + 0.000188887 dscat^3 
     * 23.59    30
     * 20.30    25
     * 16.40    20
     * 12.70    15
     * 8.60	    10
     * 6.08	    7
     * 4.36	    5
     * 2.63	    3
     * 1.75	    2
     * 1.32	    1.5
     * 0.88	    1
     */
    
    /**
     * TABLE of effective scattering thickness for 500x500 waveguides with kh=1 and holscat=dscat/500:
     * dscat    dscat_eff_approx = 1.12563 dscat + 0.00145527 dscat^2 + 0.000108717 dscat^3
     * 24.46    30
     * 16.93    20
     * 12.90    15
     * 8.72     10
     * 6.15     7
     * 4.41     5
     * 2.65     3
     * 1.77     2
     * 1.33     1.5
     * 0.89     1
     */
    
    /**
     * TABLE of effective scattering thickness for 500x500 waveguides with kh=1.05 and holscat=dscat/500:
     * dscat    dscat_eff_approx
     * 24.10	30
     * 16.66	20
     * 12.66	15
     * 8.52     10
     * 4.28     5
     */
    
    /**
     * TABLE of effective scattering thickness for 300x300 waveguides with kh=1 and holscat=dscat/300:
     * dscat    dscat_eff
     * 30       41.5617
     * 25       32.818
     * 20       25.1549
     * 15       18.246
     * 12.5     14.9555
     * 10       11.755
     * 5        5.75
     * 
     * LOCAL LINEAR FITTING:
     * .....    100
     * 34.5     50
     * 23.2     30
     * 16.27    20
     * 12.53    15
     */
    const double dscat = 0.;  // Scattering depth, L/lscat.
    const double dabso = 0.;  // Absorption depth, L/labso.
    
    const double kh = 1.;  // Wavenumber times the lattice step. Avoid kh=1 because creates resonances when ninput = 3*integer + 2.
    const double holscat = dscat/1000; // Value of h/lscat.
    const double holabso = dabso/1000; // Value of h/labso.
    const double density = 1.;  // Density of scatterers per pixel, between 0 and 1. Recommended is 1.
    
    const std::string sysname = "constriction-1c4_1000x1000/kh_" + to_string_prec(kh, 6) + "/dscat_" + to_string_prec(dscat, 6);
    
    WaveSystem sys(sysname, mesh, kh, density, holscat, holabso);
    
    //RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    //trange(0, 0) = 0.998;  trange(0, 1) = 0.002;
    //trange(1, 0) = 0.500;  trange(1, 1) = 0.01;
    //trange(2, 0) = 0.100;  trange(2, 1) = 0.005;
    //trange(3, 0) = 0.001;  trange(3, 1) = 0.00005;
    
    RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.98;  trange(0, 1) = 0.02;
    trange(1, 0) = 0.94;  trange(1, 1) = 0.02;
    trange(2, 0) = 0.90;  trange(2, 1) = 0.02;
    trange(3, 0) = 0.86;  trange(3, 1) = 0.02;
    
    //RealMatrix trange(5, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    //trange(0, 0) = 0.50;  trange(0, 1) = 0.025;
    //trange(1, 0) = 0.45;  trange(1, 1) = 0.025;
    //trange(2, 0) = 0.40;  trange(2, 1) = 0.025;
    //trange(3, 0) = 0.35;  trange(3, 1) = 0.025;
    //trange(4, 0) = 0.30;  trange(4, 1) = 0.025;
    
    //RealMatrix trange(5, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    //trange(0, 0) = 0.06;  trange(0, 1) = 0.005;
    //trange(1, 0) = 0.05;  trange(1, 1) = 0.005;
    //trange(2, 0) = 0.04;  trange(2, 1) = 0.005;
    //trange(3, 0) = 0.03;  trange(3, 1) = 0.005;
    //trange(4, 0) = 0.02;  trange(4, 1) = 0.005;
    
    Context ctx = {sys, trange};
    
    ctx.sys.printSummary();
    ctx.sys.infoHamiltonian();
    
    // Checking operations:
    //ctx.sys.setDisorder(1);
    //ctx.sys.plotMesh();
    //ctx.sys.plotMatrixHamiltonian();
    //ctx.sys.plotMatrixInputState();
    //ctx.sys.plotMatrixOutputState();
    //ctx.sys.plotGreenFunction();
    ctx.sys.plotInputModes(20);
    ctx.sys.plotTransmissionStates(20);
    //ctx.sys.checkResidual();
    //ctx.sys.checkUnitarity(true);
    
    //const int nseed = 1;    // Number of random realizations of the disorder used for averaging. Recommended for high quality: 10^4.
    //const int seed0 = 1;     // First seed used to generate realizations of the disorder. Actual seed = [seed0, seed0 + 1, ..., seed0 + nseed - 1]. 
    //                         // NB: Avoid seed0=0 for safety (some random generators are singular for seed=0).
    //const int nthread = 1;  // Number of threads used in multithreading with OpenMP.
    //const int imode = 0;     // Index of the mode in taskIMode*() and taskIAllOMP().
    
    //taskCheckUnitarityOMP(ctx.sys, nseed, seed0, nthread);
    //taskITransmissionSerial(ctx.sys, ctx.trange, nseed, seed0);
    //taskITransmissionOMP(ctx.sys, ctx.trange, nseed, seed0, nthread);
    //taskIIsotropicSerial(ctx.sys, nseed, seed0);
    //taskIIsotropicOMP(ctx.sys, nseed, seed0, nthread);
    //taskIModeOMP(ctx.sys, imode, nseed, seed0, nthread);
    //taskIAllOMP(ctx.sys, ctx.trange, imode, nseed, seed0, nthread);
    
    return 0;
}
