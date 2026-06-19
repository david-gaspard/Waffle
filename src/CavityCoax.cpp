/****
 * @date Created on 2026-06-05 at 14:53:22 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing implementation of the cavity coax object which represents a rectangular cavity coupled to a given number of coax cables.
 ***/
#include "CavityCoax.hpp"
#include "Vector2D.hpp"
#include "BaseTools.hpp"

/**
 * Constructor of the bare cavity object.
 */
CavityCoax::CavityCoax(const CavityParameters& param, const double length, const double width, 
            const int ncoaxin,  const double apercoaxin,  const double reflcoaxin, 
            const int ncoaxout, const double apercoaxout, const double reflcoaxout, 
            const int nscat, const double rscat, const double freqabso, const double h) :
        Cavity(param),
        length(length),
        width(width),
        ncoaxin(ncoaxin),
        apercoaxin(apercoaxin),
        reflcoaxin(reflcoaxin),
        ncoaxout(ncoaxout),
        apercoaxout(apercoaxout),
        reflcoaxout(reflcoaxout),
        nscat(nscat),
        rscat(rscat),
        freqabso(freqabso),
        h(h)
    {
    
    // 1. Check for possible errors:
    if (length <= 0. || width <= 0.) {
        throw std::invalid_argument("In CavityCoax(): Invalid cavity size, cannot be negative or zero.");
    }
    else if (h <= 0.) {
        throw std::invalid_argument("In CavityCoax(): Invalid mesh step, cannot be negative or zero.");
    }
    else if (ncoaxin <= 0 || ncoaxout <= 0) {
        throw std::invalid_argument("In CavityCoax(): Invalid number of coax cables, cannot be negative or zero.");
    }
    else if (apercoaxin*width <= h || apercoaxin*ncoaxin > 1.) {
        std::string msg = "In CavityCoax(): Invalid aperture of input coax cables, received apercoaxin=" + to_string_prec(apercoaxin, 6) 
                        + ", expected in [" + to_string_prec(h/width, 6) + ", " + to_string_prec(1./ncoaxin, 6) + "].";
        throw std::invalid_argument(msg);
    }
    else if (apercoaxout*width <= h || apercoaxout*ncoaxout > 1.) {
        std::string msg = "In CavityCoax(): Invalid aperture of output coax cables, received apercoaxout=" + to_string_prec(apercoaxout, 6) 
                        + ", expected in [" + to_string_prec(h/width, 6) + ", " + to_string_prec(1./ncoaxout, 6) + "].";
        throw std::invalid_argument(msg);
    }
    else if (reflcoaxin < 0. || reflcoaxin > 1. || reflcoaxout < 0. || reflcoaxout > 1.) {
        throw std::invalid_argument("In CavityCoax(): Invalid reflection probability, expected in [0, 1].");
    }
    else if (nscat < 0) {
        throw std::invalid_argument("In CavityCoax(): Invalid number of scatterers, cannot be negative.");
    }
    else if (rscat < 0.) {
        throw std::invalid_argument("In CavityCoax(): Invalid radius of scatterers, cannot be negative.");
    }
    
    // 2. Initialize the internal parameters:
    nx = static_cast<int>(std::round(length/h));  // Number of mesh point in the 'x' direction.
    ny = static_cast<int>(std::round(width/h));   // Number of mesh point in the 'y' direction.
    rscatoh = rscat/h;               // Radius of the scatterers divided by the mesh step, R_scat/h.
    holabso = 4.*PI*freqabso*h/C0;   // Value of h/l_abso, absorption per pixel.
    wcoaxin  = static_cast<int>(std::round(apercoaxin*ny)); // Number of pixels in input coax cables.
    wcoaxout = static_cast<int>(std::round(apercoaxout*ny)); // Number of pixels in output coax cables.
    
}

/**
 * Returns a summary of the present CavityCoax object.
 */
std::vector<std::string> CavityCoax::summary() const {
    std::vector<std::string> smr;
    smr.push_back("CavityCoax: Microwave rectangular cavity connected to coaxial cables and filled with metallic cylinders.");
    smr.push_back("Geometry: L=" + to_string_prec(length, 6) + "m, W=" + to_string_prec(width, 6) + "m, h=" + to_string_prec(h, 6) 
        + "m, Nx=" + std::to_string(nx) + ", Ny=" + std::to_string(ny));
    smr.push_back("Input coaxes: Ncoaxin=" + std::to_string(ncoaxin) + ", Apercoaxin=" + to_string_prec(apercoaxin, 6) 
        + ", Wcoaxin=" + std::to_string(wcoaxin) + "px, Reflcoaxin=" + to_string_prec(reflcoaxin, 6));
    smr.push_back("Output coaxes: Ncoaxout=" + std::to_string(ncoaxout) + ", Apercoaxout=" + to_string_prec(apercoaxout, 6) 
        + ", Wcoaxout=" + std::to_string(wcoaxout) + "px, Reflcoaxout=" + to_string_prec(reflcoaxout, 6));
    smr.push_back(
        "Absorption: Fabso=" + to_string_prec(freqabso, 6) + "Hz, h/labso=" + to_string_prec(holabso, 6) + ", L/labso=" + to_string_prec(nx*holabso, 6));
    
    if (getNFreq() == 1) {// When only one frequency is required, it is the center of the interval.
        const double freq  = (getFreqMin() + getFreqMax())/2;
        const int npropin  = ncoaxin  * static_cast<int>(std::floor(2.*freq*width*apercoaxin/C0));  // Number of input modes.
        const int npropout = ncoaxout * static_cast<int>(std::floor(2.*freq*width*apercoaxout/C0)); // Number of output modes.
        const double kh = 2.*PI*freq*h/C0;
        smr.push_back(
            "Frequencies: Nfreq=" + std::to_string(getNFreq()) +", Freq=" + to_string_prec(freq, 6) + "Hz, Npropin=" + std::to_string(npropin) 
            + ", Npropout=" + std::to_string(npropout) + ", kh=" + to_string_prec(kh, 6));
    }
    else {
        const int npropin_min  = ncoaxin  * static_cast<int>(std::floor(2.*getFreqMin()*width*apercoaxin/C0));  // Min number of input modes.
        const int npropin_max  = ncoaxin  * static_cast<int>(std::floor(2.*getFreqMax()*width*apercoaxin/C0));  // Max number of input modes.
        const int npropout_min = ncoaxout * static_cast<int>(std::floor(2.*getFreqMin()*width*apercoaxout/C0)); // Min number of output modes.
        const int npropout_max = ncoaxout * static_cast<int>(std::floor(2.*getFreqMax()*width*apercoaxout/C0)); // Max number of output modes.
        const double khmin = 2.*PI*getFreqMin()*h/C0;
        const double khmax = 2.*PI*getFreqMax()*h/C0;
        smr.push_back(
            "Frequencies: Nfreq=" + std::to_string(getNFreq()) +", Freq=[" + to_string_prec(getFreqMin(), 6) + "Hz, " + to_string_prec(getFreqMax(), 6) 
            + "Hz], Npropin=[" + std::to_string(npropin_min) + ", " + std::to_string(npropin_max) 
            + "], Npropout=[" + std::to_string(npropout_min) + ", " + std::to_string(npropout_max) 
            + "], kh=[" + to_string_prec(khmin, 6) + ", " + to_string_prec(khmax, 6) + "]");
    }
    smr.push_back(
        "Disorder: Nscat=" + std::to_string(nscat) + ", Rscat=" + to_string_prec(rscat, 6) + "m, Rscat/h=" + to_string_prec(rscatoh, 6) 
        + " " + (rscatoh <= 0.71 ? "(point)": "(disk)") + ", Nseed=" + std::to_string(getNSeed()) + ", Seed0=" + std::to_string(getSeed(0)));
    return smr;
}

/**
 * Generate the wave system using the given simulation parameters prescribed by the present CavityCoax object for the given "seed" used to generate
 * the realization of the disorder (random configuration of the cylinders) and the given "freq", which is the frequency (in Hz).
 */
WaveSystem CavityCoax::generateSystem(const uint64_t seed, const double freq) const {
    
    // 1. Prepare the mesh:
    SquareMesh mesh;
    mesh.addRectangle(1, nx, 1, ny, BND_MIRROR); // Add the main rectangle.
    
    // 2. Add boundary conditions to simulate coax cables:
    int y;
    if (ncoaxin == 1) {// Exception for a single input coax.
        y = static_cast<int>(std::round(1. + (ny - wcoaxin)/2.));
        mesh.setBoundaryRectangle(1, 1, y, y + wcoaxin - 1, DIR_WEST, BND_INPUT); // Setup the boundary conditions.
    }
    else {
        for (int ic = 0; ic < ncoaxin; ic++) {// Loop over input coax cables.
            y = static_cast<int>(std::round(1. + ic * static_cast<double>(ny - wcoaxin)/(ncoaxin - 1.)));
            mesh.setBoundaryRectangle(1, 1, y, y + wcoaxin - 1, DIR_WEST, BND_INPUT); // Setup the boundary conditions.
        }
    }
    if (ncoaxout == 1) {// Exception for a single output coax.
        y = static_cast<int>(std::round(1. + (ny - wcoaxout)/2.));
        mesh.setBoundaryRectangle(nx, nx, y, y + wcoaxout - 1, DIR_EAST, BND_OUTPUT); // Setup the boundary conditions.
    }
    else {
        for (int ic = 0; ic < ncoaxout; ic++) {// Loop over output coax cables.
            y = static_cast<int>(std::round(1. + ic * static_cast<double>(ny - wcoaxout)/(ncoaxout - 1.)));
            mesh.setBoundaryRectangle(nx, nx, y, y + wcoaxout - 1, DIR_EAST, BND_OUTPUT); // Setup the boundary conditions.
        }
    }
    
    // 3. Compute the random configurations of cylinders:
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
    
    // 4. Build the wave system (including the Hamiltonian):
    const double kh = 2.*PI*freq*h/C0;  // Compute the wavenumber times the mesh step from the fact that: nprop = 2*freq*width/c0 = k*width/pi = kh*ny/pi.
    const double holscat = 0.; // Scattering strength. This value is unused since setDisorder() is never called.
    const double density = 1.; // Density of the disorder (unused value).
    const double dabso = holabso * nx;  // Full absorption thickness, L/l_abso.
    const std::string sysname = "cavity-coax-in" + std::to_string(ncoaxin) + "-out" + std::to_string(ncoaxout) 
        + "_" + std::to_string(nx) + "x" + std::to_string(ny)
        + ( wcoaxin == wcoaxout ? "/wcoax_" + std::to_string(wcoaxin) : "wcoax_in_" + std::to_string(wcoaxin) + "_out_"+ std::to_string(wcoaxout) )
        + "/nscat_" + std::to_string(nscat) + "/rscatoh_" + to_string_prec(rscatoh, 6) + "/dabso_" + to_string_prec(dabso, 6);
    
    WaveSystem sys(sysname, mesh, kh, density, holscat, holabso);
    
    // 5. Add reflective barrier in front of coaxes:
    const double d2pkh2_in = laplacianEigenvalue(0, wcoaxin) + kh*kh;  // Eigenvalue of D_y^2 + (kh)^2 for the fundamental mode of the input coax.
    const dcomplex uh2_in = 2. * std::sqrt( (reflcoaxin/(1.-reflcoaxin)) * d2pkh2_in * ( 1. - d2pkh2_in/4. ) );  // Value U(x)*h^2 of the input potential barrier.
    const double d2pkh2_out = laplacianEigenvalue(0, wcoaxout) + kh*kh;  // Eigenvalue of D_y^2 + (kh)^2 for the fundamental mode of the output coax.
    const dcomplex uh2_out = 2. * std::sqrt( (reflcoaxout/(1.-reflcoaxout)) * d2pkh2_out * ( 1. - d2pkh2_out/4. ) );  // Value U(x)*h^2 of the output otential barrier.
    
    for (int y = 1; y <= ny; y++) {// Loop on the points to add to the potential barrier.
        sys.addPotential(1, y, uh2_in);
        sys.addPotential(nx, y, uh2_out);
    }
    
    return sys;
}
