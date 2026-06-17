/****
 * @date Created on 2026-05-26 at 15:28:46 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing functions for simulating a two-dimensional microwave cavity filled with cylindrical scatterers such as the one considered at ENS Rennes (Matthieu Davy's group).
 ***/
#include "CavityBare.hpp"
#include "CavityCoax.hpp"
#include "BaseTools.hpp"
#include <numeric>
#include <iomanip>

/**
 * Test the cavity for a given realization of the disorder (seed) and a given number of proapgating modes (nprop).
 */
void testCavity(const Cavity& cavity) {
    
    const uint64_t seed = 1; // Seed used for disorder generation.
    const double freq = 13.5e9; // Default frequency (Hz).
    
    WaveSystem sys = cavity.generateSystem(seed, freq);
    
    sys.printSummary(); // Print the summary of the wave system.
    //sys.plotMesh(); // Plot the mesh.
    //sys.plotMatrixHamiltonian(); // Plot the Hamiltonian matrix as a raster image.
    //sys.plotMatrixInputState(); // Plot the input state matrix as a raster image.
    //sys.plotMatrixOutputState(); // Plot the output state matrix as a raster image.
    //sys.plotGreenFunction(); // Plot all the modal components of the Green function.
    //sys.plotInputModes(3); // Plot the lowest waveguide modes.
    //sys.plotTransmissionStates(22); // Plot the lowest transmission eigenstates.
    //sys.checkResidual(); // Check the residual of the solver.
    //sys.checkUnitarity(true); // Check the unitarity of the propagation (makes sense only for no absorption).
    
    //cavity.generateSystem(2, freq).plotTransmissionStates(22);
    //cavity.generateSystem(3, freq).plotTransmissionStates(22);
    
}

/**
 * Save the data contained in "tval" (typically transmission eigenvalues) and plot the histogram.
 */
void plotHistogram(const std::vector<double>& tval, const Cavity& cavity, const int nthread, const double ctime) {
    
    const uint ntval = tval.size();  // Number of samples.
    WaveSystem sys = cavity.generateSystem(cavity.getSeed(0), cavity.getMeanFreq());
    const int prec = 16;     // Precision used in printing double precision values.
    const std::string filename = sys.uniqueFile("tspectrum", ".csv");
    const double fsize = (prec+4.) * ntval;  // Roughly estimated file size in bytes (octets).
    std::cout << TAG_INFO << "Save samples to file '" << filename << "', size ~" << (fsize/1e6) << " Mo...\n";
    
    std::ofstream ofs(filename);    // Open the output file.
    ofs << std::setprecision(prec); // Set the printing precision.
    writeTimestamp(ofs, "%% ");     // Apply a timestamp at the beginning.
    
    for (const std::string& line : cavity.summary()) {// Loop over the lines of the summary.
        ofs << "%% " << line << "\n";
    }
    const double tavg = std::reduce(tval.begin(), tval.end())/static_cast<double>(ntval);
    const double dscat_eff = 2*EXTRAPOLEN_APPROX * (1./tavg - 1.);   // Deduce the effective scattering depth (as if in a waveguide).
    ofs << "%% Results: Nsample=" << ntval << ", Tavg=" << tavg << ", L/lscat_eff=" << dscat_eff << "\n"
        << "%% Computation: Time=" << ctime << " s, Nthread=" << nthread << ".\n";
    ofs << "T\n";
    
    for (uint ival = 0; ival < ntval; ival++) {// Loop over the rows of tvalstore.
        ofs << tval.at(ival) << "\n";
    }
    ofs.close();  // Close the stream before calling an external script (this may cause I/O trouble).
    
    std::string cmd("plot/plot_histo.py " + filename);
    std::cout << TAG_EXEC << cmd << "\n";
    if (std::system(cmd.c_str())) {
        std::cout << TAG_WARN << "The plot script returned an error.\n";
    }
}

/**
 * Compute the histogram of transmission eigenvalues.
 */
void taskTSpectrum(const Cavity& cavity, const int nthread) {
    
    const int njob = cavity.getNSeed() * cavity.getNFreq();  // Total number of jobs to do.
    std::vector<double> tval_store;  // List of all transmission eigenvalues.
    // Note that transmission eigenvalues must be added dynamically because we do not know the total number in advance since the frequency is varying.
    
    int cjob = 0; // Current number of completed jobs (i.e., realizations of the disorder).
    std::string info = "TSpectrum, " + std::to_string(nthread) + " thr";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    #pragma omp parallel num_threads(nthread)
    {
        int iseed, ifreq, ntval;   // Current seed, frequency index, and number of transmission eigenvalues (OMP-private).
        
        #pragma omp for schedule(dynamic)
        for (int job = 0; job < njob; job++) {// Loop on the jobs to perform.
            
            iseed = job % cavity.getNSeed();  // Seed index (realization of the disorder).
            ifreq = job / cavity.getNSeed();  // Frequency index.
            
            WaveSystem sys = cavity.generateSystem(cavity.getSeed(iseed), cavity.getFreq(ifreq)); // Generate the wave system to solve (mesh and Hamiltonian).
            
            ntval = std::min(sys.getNInputProp(), sys.getNOutputProp());  // Number of transmission eigenvalues.
            RealMatrix tval_loc(ntval, 1);  // Allocate space for the transmission eigenvalues locally (OMP-private) because it depends on the number of modes.
            
            //std::cout << TAG_INFO << "job=" << job << ", iseed=" << iseed << ", ifreq=" << ifreq << ", ntval=" << ntval << "\n";
            
            sys.addTSpectrum(tval_loc); // Compute the transmission eigenvalues.
            
            // Critical section to gather the eigenvalues and print the progress bar:
            #pragma omp critical
            {
                for (int ival = 0; ival < ntval; ival++) {// Loop on the eigenvalues.
                    tval_store.push_back(tval_loc(ival, 0));
                }
                cjob++;
                printProgressBar(cjob, njob, info, start);
            }
        }
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    plotHistogram(tval_store, cavity, nthread, ctime);
}

/**
 * Main function of the program to simulate a two-dimensional microwave cavity filled with cylindrical scatterers.
 */
int main(int argc, char** argv) {
    
    std::cout << "****** This is " << PROGRAM_COPYRIGHT << " ******\n";
    
    // Parameters:
    const double length      = 0.400;       // Length of the cavity (in meters). Default=0.400 m.
    const double width       = 0.252;       // With of the cavity (in meters). Default=0.252 m.
    const int nscat          = 0;           // Number of scatterers. Typically: 10, 20, 50, 100, and so on.
    const double rscat       = 0.;          // Radius of the metallic scatterers (in meters). Default=0.0031 m. Zero gives single-pixel scatterers.
    const double dscat       = 0.5;         // Scattering depth, L/l_scat, of continuous Dirac-delta disorder (if present). Zero to disable.
    const double freqabso    = 0.;          // Absorption frequency (in Hz). Default=1.5e6 Hz for rscat=3.1e-3, and 7.0e6 Hz for rscat=1e-3.
    
    const int ncoaxin        = 8;           // Number of input coaxes.
    const double apercoaxin  = 1./15;       // Aperture of input coaxes (in fraction of the width). Typically, 1/15 when ncoaxin=8, and 3/8 when coaxin=1.
    
    //const int ncoaxin        = 1;           // Number of input coaxes.
    //const double apercoaxin  = 3./8;      // Aperture of input coaxes (in fraction of the width). Typically, 1/15 when ncoaxin=8, and 3/8 when coaxin=1.
    
    const int ncoaxout       = ncoaxin;     // Number of output coax cables.
    const double apercoaxout = apercoaxin;  // Aperture of output coax cables (in fraction of the width).
    
    const double h           = 5.e-4;       // Size of the step used in the space discretization (in meters). Default=5.0e-4 m.
    const int nthread        = 10;          // Number of independent computation threads used in the parallelization.
    
    const CavityParameters param = {
        .nseed = 5000,       // Number of realizations of the disorder.
        .nfreq = 1,          // Number of frequencies in the frequency range. Nfreq=1 selects the center of the interval.
        .freqmin = 12.0e9,   // Minimum field frequency of the frequency range (in Hz). Default=12.0e9 Hz.
        .freqmax = 15.0e9    // Minimum field frequency of the frequency range (in Hz). Default=15.0e9 Hz.
    };
    
    CavityBare cavity(param, length, width, nscat, rscat, dscat, freqabso, h); // Create a bare cavity with circular scatterers.
    //CavityCoax cavity(param, length, width, ncoaxin, apercoaxin, ncoaxout, apercoaxout, nscat, rscat, freqabso, h); // Create a cavity with coax guides.
    
    cavity.printSummary(); // Print the summary of the cavity parameters.
    
    testCavity(cavity);  // Test the cavity.
    
    // Compute the transmission eigenvalue spectrum:
    //taskTSpectrum(cavity, nthread);
    
    return 0;
}
