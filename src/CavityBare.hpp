/****
 * @date Created on 2026-06-04 at 17:59:48 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing header for the "CavityBare" class, a cavity which has input and output leads of same width as the actual cavity width.
 * In this regard, the cavity is more like a waveguide with one input and one output.
 ***/
#ifndef _CAVITY_BARE_H
#define _CAVITY_BARE_H
#include "Cavity.hpp"

/**
 * Abstract class defining a generic "Cavity" object.
 */
class CavityBare : public Cavity {
    
    private:
    
    double length;   // Length of the cavity (in meters). Default=0.400 m.
    double width;    // With of the cavity (in meters). Default=0.252 m.
    int nscat;       // Number of scatterers. Typically: 10, 20, 50, 100.
    double rscat;    // Radius of the metallic scatterers (in meters). Default=0.0031 m. Zero gives single-pixel scatterers.
    double freqabso; // Absorption frequency (in Hz). Default=1.5e6 Hz.
    double h;        // Size of the step used in the space discretization (in meters). Default=5.0e-4 m.
    
    // Internal parameters:
    int nx;           // Number of points of the cavity in the 'x' direction.
    int ny;           // Number of points of the cavity in the 'y' direction.
    double rscatoh;   // Radius of the scatterers divided by the mesh step, R_scat/h.
    double holabso;   // Absorption thickness per mesh point, h/l_abso.
    
    public:
    
    // Constructor:
    CavityBare(const CavityParameters& param, const double length, const double width, const int nscat, const double rscat, const double freqabso, const double h);
    
    // Methods shared by any cavity:
    std::vector<std::string> summary() const;
    WaveSystem generateSystem(const uint64_t seed, const double freq) const;
    
};

#endif
