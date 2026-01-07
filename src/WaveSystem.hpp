/****
 * @date Created on 2025-07-31 at 14:00:37 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the WaveSystem object, which describes a stationary scalar wave field defined on a SquareMesh.
 ***/
#ifndef _WAVE_SYSTEM_H
#define _WAVE_SYSTEM_H
#include "SquareMesh.hpp"
#include "SparseComplexMatrix.hpp"
#include "RealMatrix.hpp"

/**
 * Class defining the Square Mesh object.
 */
class WaveSystem {
    
    private:
    
    std::string sysname;  // String containing the name of the system (typically describing the system geometry) which is used to generate file output.
    
    SquareMesh mesh;   // Square mesh object.
    
    double kh;         // Wavenumber times the lattice step, 2*pi*h/lambda. Also the phase accumulated across a lattice step (in radian).
    double density;    // Density of scattering pixels per unit pixel. Value between 0 and 1.
    double holscat;    // Lattice step divided by the scattering mean free path, h/lscat. Total length is not a well defined unit.
    double holabso;    // Lattice step divided by the absorption length, h/labso. Total length is not a well defined unit.
    
    // TODO: Maybe store the complex wavenumber "khc" defined by kh + I*(h/labso)/2 ??
    
    int npoint;        // Total number of points in the "mesh".
    int ninput;        // Number of points in the input, also equal to the number of input modes.
    int noutput;       // Number of points in the output, also equal to the number of output modes.
    int ninputprop;    // Number of propagating modes in the input lead(s).
    int noutputprop;   // Number of propagating modes in the output lead(s).
    double dosinput;   // Density of states in the input lead(s).
    double dosoutput;  // Density of states in the input lead(s).
    double doslattice; // Free density of states on a square lattice. Exact value based on the elliptic integral K(k).
    
    bool computed;     // Flag indicating if the Green function have been computed. "true" if "green" has been computed, "false" otherwise.
                       // This flag is set to "false" each time the Hamiltonian is modified by setDisorder(), and set to "true" by computeGreenFunction().
    
    // Internal matrices:
    SparseComplexMatrix hamiltonian;  // Hamiltonian matrix, i.e., discretization of the operator (d_x^2 + d_y^2 + k^2 + i*eps - U(x,y))*h^2 
                                      // involved in the calculation of the Green function. Size: (npoint, npoint).
    SparseComplexMatrix inputState;   // Matrix storing the input modes in columns, and used as the independent term of the linear system
                                      // for the Green function. Size: (npoint, ninput).
    SparseComplexMatrix outputState;  // Matrix storing the output modes in columns. Size: (npoint, noutput).
    
    ComplexMatrix inputKlh;   // Column matrix containing the longitudinal wavenumbers in each channel. Size: (ninput, 1).
    ComplexMatrix outputKlh;  // Column matrix containing the longitudinal wavenumbers in each channel. Size: (noutput, 1).
    ComplexMatrix green;  // Matrix containing in columns the Green functions associated to each input modes. Size: (npoint, ninput).
    
    public:
    
    // Constructors/Destructors:
    WaveSystem(const std::string& name, const SquareMesh& mesh, const double kh, const double density, const double holscat, const double holabso);
    
    // Getters:
    int getNPoint() const;
    MeshPoint getPoint(const unsigned int ipoint) const;
    int getNInput() const;
    int getNOutput() const;
    int getNInputProp() const;
    int getNOutputProp() const;
    double getWavenumber() const;
    double getDensity() const;
    double getScattering() const;
    double getAbsorption() const;
    std::string getName() const;
    std::string uniqueFile(const std::string& dataname, const std::string& extension) const;
    std::vector<std::string> summary() const;
    
    // Print/save methods:
    void printSummary() const;
    void infoHamiltonian() const;
    void plotHamiltonian() const;
    void plotInputState() const;
    void plotOutputState() const;
    void plotMesh() const;
    void plotIntensity(const RealMatrix& intensity, const std::string& description, const std::string& filename) const;
    void plotGreenFunction();
    void plotTransmissionStates();
    
    // Public computational methods:
    void setDisorder(const uint64_t seed);
    void transmissionMatrix(ComplexMatrix& tmat);
    void reflectionMatrix(ComplexMatrix& rmat);
    void checkUnitarity(const bool showtval);
    void checkResidual();
    void addITransmission(const RealMatrix& trange, RealMatrix& tprofile, RealMatrix& nsample, RealMatrix& tval);
    void addIIsotropic(RealMatrix& iavg);
    void addIMode(RealMatrix& iavg, const int imode);
    
    private:
    
    // Private setters (because the Hamiltonian should be recomputed afterwards in principle):
    double checkWavenumber(const double kh) const;
    double checkDensity(const double density) const;
    double checkScattering(const double holscat) const;
    double checkAbsorption(const double holabso) const;
    
    // Private computational methods:
    void computeHamiltonian();
    int computeNInputProp() const;
    int computeNOutputProp() const;
    void computeIOStates();
    void computeGreenFunction();
    void computeTransmissionStates(ComplexMatrix& tstate, RealMatrix& tval);
    
};

#endif
