/****
 * @date Created on 2026-05-26 at 15:28:46 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing functions for simulating a two-dimensional microwave cavity filled with cylindrical scatterers such as the one considered at ENS Rennes (Matthieu Davy's group).
 ***/
#include "WaveSystem.hpp"
#include "BaseTools.hpp"
#include <random>
#include <iomanip>

/**
 * Defines all the simulation parameters for the bare cavity.
 * Note that these parameters must be dimensionless (no dimensionfull parameter allowed).
 */
struct CavityBare {
    
    std::string name; // Name of the cavity.
    int nx;           // Number of points of the cavity in the 'x' direction.
    int ny;           // Number of points of the cavity in the 'y' direction.
    int nscat;        // Number of metallic cylindrical scatterers.
    double rscatoh;   // Radius of the scatterers divided by the mesh step, R_scat/h.
    double holabso;   // Absorption thickness per mesh point, h/l_abso.
    double npropfmin; // Minimum number of propagating modes from the frequency range (npropfmin = kmin*width/pi).
    double npropfmax; // Maximum number of propagating modes from the frequency range (npropfmax = kmax*width/pi).
    int nfreq;        // Number of requencies over which to average. The interval [npropfmin, npropfmax] is covered uniformly at constant steps.
    int nseed;        // Number of realizations of the disorder.
    int seed0;        // Initial seed of the disorder generation.
    
    std::vector<std::string> summary() const;
    void printSummary() const;
};

/**
 * Returns a summary.
 */
std::vector<std::string> CavityBare::summary() const {
    std::vector<std::string> smr;
    smr.push_back("Transmission eigenvalues in a microwave cavity with Nscat metallic cylinders of radius Rscat.");
    smr.push_back("Bare cavity: Nx=" + std::to_string(nx) + ", Ny=" + std::to_string(ny) + ", h/labso=" + std::to_string(holabso));
    smr.push_back("Frequencies: Nfreq=" + std::to_string(nfreq) +" in Nprop=kW/pi=[" + std::to_string(npropfmin) + ", " + std::to_string(npropfmax) + "]");
    smr.push_back("Disorder: Nscat=" + std::to_string(nscat) + ", Rscat/h=" + std::to_string(rscatoh) + " " + (rscatoh <= 0.71 ? "(point)": "(disk)")
        + ", Nseed=" + std::to_string(nseed) + ", Seed0=" + std::to_string(seed0));
    return smr;
}

/**
 * Prints the summary.
 */
void CavityBare::printSummary() const {
    for (const std::string& line : summary()) {// Loop over the lines of the summary.
        std::cout << TAG_INFO << line << "\n";
    }
}

/**
 * Generate a list of positions of random non-overlapping disks of given "radius" in the rectangle xmin, xmax, ymin, ymax, using the given "seed".
 * The positions avoid the edges of the rectangle.
 */
std::vector<Vector2D> nonOverlappingDisks(const int ndisk, const double radius,
        const double xmin, const double xmax, const double ymin, const double ymax, const uint64_t seed) {
    
    // 1. Check for possible errors:
    if (ndisk < 0) {
        throw std::invalid_argument("In nonOverlappingDisks(): Invalid number of disks, cannot be negative.");
    }
    else if (radius < 0.) {
        throw std::invalid_argument("In nonOverlappingDisks(): Invalid disk radius, cannot be negative.");
    }
    else if (xmax <= xmin || ymax <= ymin) {
        throw std::invalid_argument("In nonOverlappingDisks(): Invalid rectangular bounds, xmax <= xmin or ymax <= ymin.");
    }
    
    const double pfrac = PI/(6.*std::sqrt(3.));  // Minimum packing fraction leading to congestion (about 30%).
    const double surf = std::max(xmax - xmin - 2.*radius, 0.) * std::max(ymax - ymin - 2.*radius, 0.);  // Available surface.
    const int ndisk_max = static_cast<int>(std::floor(surf*pfrac/(PI*radius*radius)));
    if (ndisk > ndisk_max && ndisk_max > 0) {// Check if the number of disk is not too large, and if the maximum ndisk_max is correct.
        std::string msg = "In nonOverlappingDisks(): Number of disks (ndisk=" + std::to_string(ndisk) + ") is too large, maximum is " 
            + std::to_string(ndisk_max) + " (max packing is " + std::to_string(100.*pfrac) + "%).";
        throw std::invalid_argument(msg);
    }
    
    // 2. Prepare the random generators:
    std::mt19937_64 rng;  // Instantiate the standard Mersenne Twister random number generator (64-bit-return version).
    rng.seed(seed);       // Initialize the random generator with the given seed.
    std::uniform_real_distribution<double> random_uniform_x(xmin + radius, xmax - radius);  // Avoid the edges.
    std::uniform_real_distribution<double> random_uniform_y(ymin + radius, ymax - radius);
    
    double x, y, dx, dy;
    bool valid;
    std::vector<Vector2D> poslist;
    
    for (int i = 0; i < ndisk; i++) {// Loop over the disks to add.
        
        do {// Loop over trials.
            x = random_uniform_x(rng);  // Generate a random position uniformly in the allowed rectangle.
            y = random_uniform_y(rng);
            
            // Check if the disk does not overlap another disk:
            valid = true;
            for (uint j = 0 ; j < poslist.size(); j++) {
                dx = x - poslist.at(j).x;
                dy = y - poslist.at(j).y;
                if (dx*dx + dy*dy < 4*radius*radius) {// Overlap detected.
                    valid = false;
                    break;
                }
            }
            
        } while (not valid);
        
        poslist.push_back(Vector2D(x, y)); // Save the disk to check the overlaps.
    }
    
    return poslist;
}

/**
 * Generate the wave system using the given simulation parameters prescribed by "cavity" and the specific "seed" used to generate
 * the realization of the disorder (random configuration of the cylinders) and the specific "nprop", which is the number of propagating modes
 * given by nprop = 2*freq*width/c0 = k*width/pi = kh*ny/pi.
 */
WaveSystem generateSystem(const CavityBare& cavity, const uint64_t seed, const double nprop) {
    
    if (VERBOSE >= 1) {
        std::cout << TAG_INFO << "Generate WaveSystem for seed=" << seed << ", nprop=" << nprop << "...\n";
    }
    
    // 1. Prepare the mesh:
    SquareMesh mesh;
    mesh.addRectangle(1, cavity.nx, 1, cavity.ny, BND_MIRROR); // Add the main rectangle.
    mesh.setBoundaryRectangle(1, 1, 1, cavity.ny, DIR_WEST, BND_INPUT); // Setup the boundary conditions.
    mesh.setBoundaryRectangle(cavity.nx, cavity.nx, 1, cavity.ny, DIR_EAST, BND_OUTPUT);
    
    // 2. Compute the random configurations of cylinders:
    std::vector<Vector2D> poslist = nonOverlappingDisks(cavity.nscat, cavity.rscatoh, 0.5, cavity.nx+0.5, 0.5, cavity.ny+0.5, seed);
    for (const Vector2D& p : poslist) {// Loop on the positions.
        if (cavity.rscatoh <= 0.71) {// If the scatterers are too small, then consider that it is more appropriate to punch pixels.
            mesh.removePoint(p.x, p.y); // Remove a single pixel.
        }
        else {
            mesh.removeDisk(p.x, p.y, cavity.rscatoh); // Punch holes in the mesh to create the scatterers.
        }
    }
    mesh.finalize(); // Finalize the mesh (calculate the neighbors).
    
    // 3. Build the wave system (including the Hamiltonian):
    const double holscat = 1.; // Scattering strength. This value is unused since setDisorder() is never called.
    const double density = 1.; // Density of the disorder (unused value).
    const double kh = PI*nprop/cavity.ny;  // Compute the wavenumber times the mesh step from the fact that: nprop = 2*freq*width/c0 = k*width/pi = kh*ny/pi.
    const double dabso = cavity.holabso * cavity.nx;  // Full absorption thickness, L/l_abso.
    const std::string sysname = cavity.name + "_" + std::to_string(cavity.nx) + "x" + std::to_string(cavity.ny)
        + "/nscat_" + std::to_string(cavity.nscat) + "/rscatoh_" + to_string_prec(cavity.rscatoh, 6) + "/dabso_" + to_string_prec(dabso, 6);
    
    return WaveSystem(sysname, mesh, kh, density, holscat, cavity.holabso);
}

/**
 * Test the cavity for a given realization of the disorder (seed) and a given number of proapgating modes (nprop).
 */
void testSystem(const CavityBare& cavity, const uint64_t seed, const double nprop) {
    
    WaveSystem sys = generateSystem(cavity, seed, nprop);
    
    // Common tests:
    //sys.plotMesh(); // Plot the mesh.
    //sys.plotMatrixHamiltonian(); // Plot the Hamiltonian matrix as a raster image.
    //sys.plotMatrixInputState(); // Plot the input state matrix as a raster image.
    //sys.plotMatrixOutputState(); // Plot the output state matrix as a raster image.
    //sys.plotGreenFunction(); // Plot all the modal components of the Green function.
    //sys.plotInputModes(3); // Plot the lowest waveguide modes.
    sys.plotTransmissionStates(22); // Plot the lowest transmission eigenstates.
    //sys.checkResidual(); // Check the residual of the solver.
    //sys.checkUnitarity(true); // Check the unitarity of the propagation (makes sense only for no absorption).
}

/**
 * Save the data contained in "tval" (typically transmission eigenvalues) and plot the histogram.
 */
void plotHistogram(const std::vector<double>& tval, const CavityBare& cavity, const int nthread, const double ctime) {
    
    const uint ntval = tval.size();  // Number of samples.
    WaveSystem sys = generateSystem(cavity, cavity.seed0, (cavity.npropfmin + cavity.npropfmax)/2.);
    const char* sep = ", ";  // Separator used between entries of the CSV file.
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
    
    for (int ival = 0; ival < ntval; ival++) {// Loop over the rows of tvalstore.
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
void taskTSpectrum(const CavityBare& cavity, const int nthread) {
    
    const int njob = cavity.nseed * std::max(cavity.nfreq, 1);  // Total number of jobs to do.
    std::vector<double> tval;  // List of all transmission eigenvalues.
    
    uint64_t seed; // Current seed.
    double nprop;  // Current number of propagating modes, value of 2*freq*width/c0 = k*width/pi = kh*ny/pi.
    
    int cjob = 0; // Current number of completed jobs (i.e., realizations of the disorder).
    std::string info = "TSpectrum, " + std::to_string(nthread) + " thr";
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    #pragma omp parallel num_threads(nthread)
    {
        std::vector<double> tval_loc; // List of transmission eigenvalues (OMP-private).
        
        #pragma omp for schedule(dynamic)
        for (int job = 0; job < njob; job++) {// Loop on the jobs to perform.
            
            seed = cavity.seed0 + job % cavity.nseed;
            if (cavity.nfreq <= 1) {// If only one frequency, then take the central value.
                nprop = (cavity.npropfmin + cavity.npropfmax)/2;
            }
            else {
                nprop = cavity.npropfmin + (cavity.npropfmax - cavity.npropfmin) * static_cast<double>(job/cavity.nseed)/(cavity.nfreq - 1.);
            }
            
            WaveSystem sys = generateSystem(cavity, seed, nprop); // Generate the wave system (OMP-private).
            
            // TODO: WARNING: There are still some race conditions, use std::vector<RealMatrix> for "tval" instead and a "critical" section.
            
            sys.addTSpectrum(tval_loc); // Compute the transmission eigenvalues.
            
            // Critical section to print the progress bar:
            #pragma omp critical
            {
                cjob++;
                printProgressBar(cjob, njob, info, start);
            }
        }
        
        // Critical section to gather all data together:
        #pragma omp critical
        {
            tval.insert(tval.end(), tval_loc.begin(), tval_loc.end());
        }
    }
    
    // End progress bar and compute the total time (in seconds):
    const double ctime = endProgressBar(start);
    
    plotHistogram(tval, cavity, nthread, ctime);
}

/**
 * Main function of the program to simulate a two-dimensional microwave cavity filled with cylindrical scatterers.
 */
int main(int argc, char** argv) {
    
    std::cout << "****** This is " << PROGRAM_COPYRIGHT << " ******\n";
    
    // Parameters:
    const std::string name = "cavity-bare";  // Name of the system.
    const double c0 = 299792458.0;   // Speed of light in vacuum (299 792 458 m/s), definition of the meter.
    const double length   = 0.400;   // Length of the cavity (in meters). Default=0.400 m.
    const double width    = 0.252;   // With of the cavity (in meters). Default=0.252 m.
    const double rscat    = 3.1e-3;  // Radius of the metallic scatterers (in meters). Default=0.0031 m. Zero gives single-pixel scatterers.
    const double freqmin  = 13.5e9;  // Minimum frequency of the frequency range (in Hz). Default=12.0e9 Hz.
    const double freqmax  = 13.5e9;  // Minimum frequency of the frequency range (in Hz). Default=15.0e9 Hz.
    const double freqabso = 0.;      // Absorption frequency (in Hz). Default=1.5e6 Hz.
    const double h        = 5.e-4;   // Size of the step used in the space discretization (in meters). Default=5.0e-4 m.
    const int nscat = 10;            // Number of scatterers. Typically: 10, 20, 50, 100.
    const int nfreq = 1;             // Number of frequencies in the frequency range. Nfreq=1 selects the center of the interval.
    const int nseed = 1000;          // Number of realizations of the disorder. Recommended is 1000.
    const uint64_t seed0 = 1;        // Initial seed used in the disorder generation.
    const int nthread = 10;          // Number of independent computation threads used in the parallelization.
    
    const CavityBare cavity = {// Initialize the dimensionless simulation parameters.
        .name = name,   // Name of the system (used in file output).
        .nx = static_cast<int>(std::round(length/h)),  // Number of mesh point in the 'x' direction.
        .ny = static_cast<int>(std::round(width/h)),   // Number of mesh point in the 'y' direction.
        .nscat = nscat,                   // Number of scatterers.
        .rscatoh = rscat/h,               // Radius of the scatterers divided by the mesh step, R_scat/h.
        .holabso = 4.*PI*freqabso*h/c0,   // Value of h/l_abso, absorption per pixel.
        .npropfmin = 2*freqmin*width/c0,  // Minimum number of propagating modes from the frequency range (npmin = kmin*width/pi).
        .npropfmax = 2*freqmax*width/c0,  // Maximum number of propagating modes from the frequency range (npmin = kmin*width/pi).
        .nfreq = nfreq,  // Number of frequencies over which to average. The interval [npropfmin, npropfmax] is covered uniformly at constant steps.
        .nseed = nseed,  // Number of realizations of the disorder over which to average.
        .seed0 = seed0   // Initialize seed used by the disorder generation.
    };
    
    cavity.printSummary(); // Print the summary of the cavity parameters.
    
    // Test the system:
    const uint64_t seed = 1; // Seed used for disorder generation.
    const double nprop = 22.5; // Number of propagating modes.
    testSystem(cavity, seed, nprop);
    
    // Compute the transmission eigenvalue spectrum:
    //taskTSpectrum(cavity, nthread);
    
    return 0;
}
