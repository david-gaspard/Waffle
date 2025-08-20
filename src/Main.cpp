/****
 * @date Created on 2025-08-13 at 12:43:52 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the main functions of the Waffle program.
 ***/
#include "WaveSystem.hpp"
#include "BaseTools.hpp"

/**
 * @todo TODO:
 * DONE: (1) In WaveSystem: If possible, compute the exact expression of the free DOS on a square lattice. It is used in setDisorder()...
 * (2): In WaveSystem: Check in doing math that the normalization of transmission eigenstates is correct.
 * (3) In Main: Add histogram of transmission eigenvalues...
 * (4) In Main: Add parallelized taskTransmissionOMP()...
 * (5) In plot_map.py: Fix the partial read bug with the CSV reader which occurs when it is called by the C++ program...
 * (6) In SquareMesh: Create addImage(filename) to import PNG in order to facilitate the mesh creation...
 * (7) In SparseComplexMatrix: Improve plotImage() to export directly a PNG image instead of a heavy PPM file.
 */

/**
 * Create a waveguide-shape system of given "length" and "width" (in number of lattice points),
 * and defines the range of transmission values of interest for computing the disorder-averaged profile of transmission eigenchannels. 
 */
WaveSystem createWaveguide() {
    
    const int length = 120;   // Number of lattice points in the longitudinal direction, L/h.
    const int width  = 120;   // Number of lattice points in the transverse direction, W/h.
    
    const double dscat = 5.;  // Scattering depth, L/lscat.
    const double dabso = 0.;  // Absorption depth, L/labso.
    
    const std::string name = "Waveguide_" + std::to_string(length) + "x" + std::to_string(width) + "/dscat_" + to_string_prec(dscat, 6);
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, length, 0, width, BND_MIRROR);
    //mesh.removeDisk(length/3, width/3, 2.5);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, 0, width, DIR_WEST, BND_INPUT);
    //mesh.setBoundaryRectangle(0, 0, width/4, 3*width/4, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(length, length, 0, width, DIR_EAST, BND_OUTPUT);
    //mesh.setBoundaryRectangle(length, length, width/4, 3*width/4, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double kh = 1.;  // Wavenumber times the lattice step. Recommended: kh = 1 -> lambda/h = 6.
    const double holscat = dscat/length;
    const double holabso = dabso/length;
    
    return WaveSystem(name, mesh, kh, holscat, holabso);
}

/**
 * Create a slab system of given length and width.
 */
WaveSystem createSlab() {
    
    const int length  = 80;   // Number of lattice points in the longitudinal direction, L/h.
    const int width   = 240;  // Number of lattice points in the transverse direction, W/h.
    const int winput  = 120;  // Width of the input channel.
    const int woutput = 120;  // Width of the output channel.
    
    const double dscat = 5.;  // Scattering depth, L/lscat.
    const double dabso = 0.;  // Absorption depth, L/labso.
    
    const std::string name = "Slab_" + std::to_string(length) + "x" + std::to_string(width) + "/dscat_" + to_string_prec(dscat, 6);
    
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
    
    return WaveSystem(name, mesh, kh, holscat, holabso);
}

/**
 * Compute all quantities related to transmission for the given WaveSystem "sys", in particular the transmission-eigenvalue distribution,
 * and the intensity profile of selected transmission eigenchannels (aka eigenstates).
 * Save the results to files with automated names, and call external plot scripts.
 */
void taskTransmissionSerial(WaveSystem& sys, RealMatrix& trange) {
    
    const uint64_t nseed = 200;  // Number of realizations of the disorder.
    
    const int npoint = sys.getNPoint();
    const int nprofile = trange.getNrow();  // Get the desired number of profiles.
    const int ntval = std::min(sys.getNInput(), sys.getNOutput());  // Get the number of transmission eigenvalues.
    const int nprop = std::min(sys.getNInputProp(), sys.getNOutputProp());  // Get the number of propagationg modes (nonzero transmission eigenvalues).
    RealMatrix tprofile(npoint, nprofile), nsample(nprofile, 1), tval(ntval, 1);
    double tmax, tmax_max = 0., tmax_min = 1., tavg = 0.;
    
    std::string msg = "Tprofile, serial";
    
    auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (uint64_t seed = 1; seed <= nseed; seed++) {// Loop over realizations of the disorder.
        
        sys.setDisorder(seed);
        sys.transmissionProfiles(trange, tprofile, nsample, tval);
        
        // TODO: Create a histogram of transmission eigenvalues "tval"......
        
        tavg += tval.sum();  // Add all transmission eigenvalues.
        
        tmax = tval.max();
        if (tmax > tmax_max) tmax_max = tmax;
        if (tmax < tmax_min) tmax_min = tmax;
        
        printProgressBar(seed, nseed, msg, start);
    }
    tavg /= nprop*nseed;  // Normalize the average transmission probability.
    
    // End progress bar and compute the total time (in seconds):
    double ctime = endProgressBar(start);
    
    // Compute the average profiles by dividing by the number of samples:
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Loop over the transmission eigenstate profiles.
        if (nsample(iprofile, 0) > MEPS) {// Only normalize the profile if the number of samples is nonzero.
            for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points.
                tprofile(ipoint, iprofile) /= nsample(iprofile, 0);
            }
        }
    }
    
    // Save the transmission eigenstate profiles to a file and plot them:
    std::string info = "Transmission eigenstate profiles for Trange=[";
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Save the ranges of transmission eigenvalues for which 
        info += " " + std::to_string(trange(iprofile, 0)) + "+-" + std::to_string(trange(iprofile, 1)) + " ";
    }
    info += "], found Nsample=[";
    for (int iprofile = 0; iprofile < nprofile; iprofile++) {// Save the ranges of transmission eigenvalues for which 
        info += " " + std::to_string(nsample(iprofile, 0)) + " ";
    }
    info += "], Nseed=" + std::to_string(nseed) + ", Tavg=" + std::to_string(tavg) + ", Tmax=[" + std::to_string(tmax_min) + ", " 
          + std::to_string(tmax_max) + "], Computation_time=" + std::to_string(ctime) + " s.";
    sys.plotIntensity(tprofile, info, "tprofile");  // Plot the average transmission eigenstate profiles.
}

/**
 * Main function of the program.
 */
int main(int argc, char** argv) {
    
    std::cout << "****** This is " << PROGRAM_COPYRIGHT << " ******\n";
    
    //// Create the wave system:
    //WaveSystem sys = createWaveguide();
    //RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    //trange(0, 0) = 0.980;  trange(0, 1) = 0.020;
    //trange(1, 0) = 0.500;  trange(1, 1) = 0.050;
    //trange(2, 0) = 0.100;  trange(2, 1) = 0.020;
    //trange(3, 0) = 0.005;  trange(3, 1) = 0.001;
    
    // Create the wave system:
    WaveSystem sys = createSlab();
    RealMatrix trange(4, 2); // Defines the selected transmission intervals for computing the averaged profile of transmission eigenchannels:
    trange(0, 0) = 0.800;  trange(0, 1) = 0.050;
    trange(1, 0) = 0.300;  trange(1, 1) = 0.020;
    trange(2, 0) = 0.100;  trange(2, 1) = 0.010;
    trange(3, 0) = 0.005;  trange(3, 1) = 0.001;
    
    sys.printInfo();
    //sys.infoHamiltonian();
    //sys.plotMesh();
    //sys.plotInputState();
    //sys.plotOutputState();
    //sys.plotHamiltonian();
    //sys.plotGreenFunction();
    //sys.checkUnitarity(true);
    
    //sys.setDisorder(1);
    //sys.plotGreenFunction();
    //sys.plotTransmissionStates();
    
    //taskTransmissionSerial(sys, trange);
    
    return 0;
}
