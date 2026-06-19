/****
 * @date Created on 2026-06-05 at 14:31:36 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing header for the cavity coax object, a cavity connected to a given number of coax cables to input and output leads.
 ***/
#ifndef _CAVITY_COAX_H
#define _CAVITY_COAX_H
#include "Cavity.hpp"

/**
 * Abstract class defining a generic "Cavity" object.
 */
class CavityCoax : public Cavity {
    
    private:
    
    double length;      // Length of the cavity (in meters). Default=0.400 m.
    double width;       // With of the cavity (in meters). Default=0.252 m.
    
    int ncoaxin;        // Number of input coax cables.
    double apercoaxin;  // Aperture of input coax cables (in fraction of the width).
    double reflcoaxin;  // Reflectivity of input coax cables, probability taking into account impedance mismatch (between 0 and 1).
    
    int ncoaxout;       // Number of output coax cables.
    double apercoaxout; // Aperture of output coax cables.
    double reflcoaxout;  // Reflectivity of output coax cables, probability taking into account impedance mismatch (between 0 and 1).
    
    int nscat;          // Number of scatterers. Typically: 10, 20, 50, 100.
    double rscat;       // Radius of the metallic scatterers (in meters). Default=0.0031 m. Zero gives single-pixel scatterers.
    double freqabso;    // Absorption frequency (in Hz). Default=1.5e6 Hz.
    double h;           // Size of the step used in the space discretization (in meters). Default=5.0e-4 m.
    
    // Internal parameters:
    int nx;             // Number of points of the cavity in the 'x' direction.
    int ny;             // Number of points of the cavity in the 'y' direction.
    int wcoaxin;        // Number of pixels of input coax cables, which is fixed (roundoff are pushed in the intervals between coaxes).
    int wcoaxout;       // Number of pixels of output coax cables, which is fixed (roundoff are pushed in the intervals between coaxes).
    double rscatoh;     // Radius of the scatterers divided by the mesh step, R_scat/h.
    double holabso;     // Absorption thickness per mesh point, h/l_abso.
    
    public:
    
    // Constructor:
    CavityCoax(const CavityParameters& param, const double length, const double width, 
            const int ncoaxin,  const double apercoaxin,  const double reflcoaxin, 
            const int ncoaxout, const double apercoaxout, const double reflcoaxout, 
            const int nscat, const double rscat, const double freqabso, const double h);
    
    // Methods shared by any cavity:
    std::vector<std::string> summary() const;
    WaveSystem generateSystem(const uint64_t seed, const double freq) const;
    
};

#endif
