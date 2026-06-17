/****
 * @date Created on 2026-06-04 at 18:23:12 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing implementation of the CavityBare class.
 ***/
#include "CavityBare.hpp"
#include "Vector2D.hpp"
#include "BaseTools.hpp"

/**
 * Constructor of the bare cavity object.
 */
CavityBare::CavityBare(const CavityParameters& param, const double length, const double width, const int nscat, const double rscat, const double dscat, const double freqabso, const double h) :
        Cavity(param),
        length(length),
        width(width),
        nscat(nscat),
        rscat(rscat),
        dscat(dscat),
        freqabso(freqabso),
        h(h)
    {
    
    // 1. Check for possible errors:
    if (length <= 0. || width <= 0.) {
        throw std::invalid_argument("In CavityBare(): Invalid cavity size, cannot be negative or zero.");
    }
    else if (h <= 0.) {
        throw std::invalid_argument("In CavityBare(): Invalid mesh step, cannot be negative or zero.");
    }
    else if (nscat < 0) {
        throw std::invalid_argument("In CavityBare(): Invalid number of scatterers, cannot be negative.");
    }
    else if (rscat < 0.) {
        throw std::invalid_argument("In CavityBare(): Invalid radius of scatterers, cannot be negative.");
    }
    
    // 2. Initialize the internal parameters:
    nx = static_cast<int>(std::round(length/h));  // Number of mesh point in the 'x' direction.
    ny = static_cast<int>(std::round(width/h));   // Number of mesh point in the 'y' direction.
    rscatoh = rscat/h;               // Radius of the scatterers divided by the mesh step, R_scat/h.
    holabso = 4.*PI*freqabso*h/C0;   // Value of h/l_abso, absorption per pixel.
    
}

/**
 * Returns a summary of the present CavityBare object.
 */
std::vector<std::string> CavityBare::summary() const {
    std::vector<std::string> smr;
    smr.push_back("CavityBare: Microwave rectangular waveguide with Nscat metallic cylinders of radius Rscat.");
    smr.push_back("Geometry: L=" + to_string_prec(length, 6) + "m, W=" + to_string_prec(width, 6) + "m, h=" + to_string_prec(h, 6) 
        + "m, Nx=" + std::to_string(nx) + ", Ny=" + std::to_string(ny));
    smr.push_back("Absorption: Fabso=" + to_string_prec(freqabso, 6) + "Hz, h/labso=" + to_string_prec(holabso, 6) + ", L/labso=" + to_string_prec(nx*holabso, 6));
    if (getNFreq() == 1) {// When only one frequency is required, it is the center of the interval.
        const double freq = (getFreqMin() + getFreqMax())/2;
        const double nprop = 2.*freq*width/C0;
        const double kh = 2.*PI*freq*h/C0;
        smr.push_back("Frequencies: Nfreq=" + std::to_string(getNFreq()) +", Freq=" + to_string_prec(freq, 6) + "Hz, Nprop=kW/pi=" + to_string_prec(nprop, 6) + ", kh=" + to_string_prec(kh, 6));
    }
    else {
        const double npropmin = 2.*getFreqMin()*width/C0;
        const double npropmax = 2.*getFreqMax()*width/C0;
        const double khmin = 2.*PI*getFreqMin()*h/C0;
        const double khmax = 2.*PI*getFreqMax()*h/C0;
        smr.push_back("Frequencies: Nfreq=" + std::to_string(getNFreq()) +", Freq=[" + to_string_prec(getFreqMin(), 6) + "Hz, " + to_string_prec(getFreqMax(), 6) 
            + "Hz], Nprop=kW/pi=[" + to_string_prec(npropmin, 6) + ", " + to_string_prec(npropmax, 6) 
            + "], kh=[" + to_string_prec(khmin, 6) + ", " + to_string_prec(khmax, 6) + "]");
    }
    smr.push_back("Disorder: Nscat=" + std::to_string(nscat) + ", Rscat=" + to_string_prec(rscat, 6) + "m, Rscat/h=" + to_string_prec(rscatoh, 6) 
        + " " + (rscatoh <= 0.71 ? "(point)": "(disk)") + (dscat != 0. ? ", additional continuous disorder: L/lscat=" + to_string_prec(dscat, 6) : "")
        + ", Nseed=" + std::to_string(getNSeed()) + ", Seed0=" + std::to_string(getSeed(0)));
    return smr;
}

/**
 * Generate the wave system using the given simulation parameters prescribed by the present CavityBare object, at the given "seed" used to generate
 * the realization of the disorder (random configuration of the cylinders) and the given "freq", which is the frequency (in Hz).
 */
WaveSystem CavityBare::generateSystem(const uint64_t seed, const double freq) const {
    
    // 1. Prepare the mesh:
    SquareMesh mesh;
    mesh.addRectangle(1, nx, 1, ny, BND_MIRROR); // Add the main rectangle.
    mesh.setBoundaryRectangle(1, 1, 1, ny, DIR_WEST, BND_INPUT); // Setup the boundary conditions.
    mesh.setBoundaryRectangle(nx, nx, 1, ny, DIR_EAST, BND_OUTPUT);
    
    // 2. Compute the random configurations of cylinders:
    std::vector<Vector2D> poslist = nonOverlappingDisks(nscat, rscatoh, 1.5, nx - 0.5, 1.5, ny - 0.5, seed);
    for (const Vector2D& p : poslist) {// Loop on the positions.
        if (rscatoh <= 0.71) {// If the scatterers are too small, then consider that it is more appropriate to punch pixels.
            mesh.removePoint(static_cast<int>(std::round(p.x)), static_cast<int>(std::round(p.y))); // Remove a single pixel.
        }
        else {
            mesh.removeDisk(p.x, p.y, rscatoh); // Punch holes in the mesh to create the scatterers.
        }
    }
    mesh.finalize(); // Finalize the mesh (calculate the neighbors).
    
    // 3. Build the wave system (including the Hamiltonian):
    const double kh = 2.*PI*freq*h/C0;  // Compute the wavenumber times the mesh step from the fact that: nprop = 2*freq*width/c0 = k*width/pi = kh*ny/pi.
    const double holscat = dscat/nx; // Scattering strength.
    const double density = 1.; // Density of the disorder (unused value).
    const double dabso = holabso * nx;  // Full absorption thickness, L/l_abso.
    const std::string sysname = "cavity-bare_" + std::to_string(nx) + "x" + std::to_string(ny)
        + "/nscat_" + std::to_string(nscat) + "/rscatoh_" + to_string_prec(rscatoh, 6) + "/dabso_" + to_string_prec(dabso, 6);
    
    WaveSystem sys(sysname, mesh, kh, density, holscat, holabso);
    sys.setDisorder(seed); // Add a continuous disorder background.
    return sys;
}
