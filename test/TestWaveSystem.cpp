/****
 * @date Created on 2025-08-09 at 17:52:49 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the WaveSystem object.
 ***/
#include "WaveSystem.hpp"
#include <iostream>

/**
 * Create a test system.
 */
WaveSystem createTestSystem() {
    
    SquareMesh mesh; // Declare the mesh.
    
    // Setup the mesh from geometry:
    mesh.addRectangle(-10, 10, -5, 5, BND_MIRROR);
    mesh.addDisk(10, 5, 7.5, BND_MIRROR);
    mesh.removeDisk(10, 5, 2.5);
    
    // Setup boundary conditions (only vacuum because mirror is default):
    mesh.setBoundaryRectangle(-10, -10, -5, 5, DIR_WEST, BND_OPEN);
    mesh.setBoundaryDisk(-10, 0, 2.5, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(-3, 5, 4.5, DIR_NORTH, BND_INPUT);
    mesh.setBoundaryRectangle( 10,  10, -5, 5, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryDisk(10, 5, 3.5, DIR_ALL, BND_OPEN);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.finalize();
    
    double kh = PI/2;
    double holscat = 0.1;
    double holabso = 0.;
    double density = 1.;
    
    WaveSystem sys("TestHamiltonian", mesh, kh, density, holscat, holabso);
    
    return sys;
}

/**
 * Test the construction of the Hamiltonian.
 */
int testHamiltonian() {
    std::cout << "====== TEST HAMILTONIAN AND INDEPENDENT TERMS ======\n";
    
    WaveSystem sys = createTestSystem();
    
    sys.infoHamiltonian();  // Summary of the Hamiltonian to stdout.
    sys.plotMatrixHamiltonian();  // Plot the Hamiltonian matrix to an image file.
    sys.plotMatrixInputState();   // Plot the input state matrix to an image file.
    sys.plotMatrixOutputState();  // Plot the output state matrix to a file.
    sys.plotMesh();         // Plot the mesh for checking.
    
    return 0;
}

/**
 * Create a test system for unitarity (probability conservation) checking.
 */
WaveSystem createUnitarySystem() {
    
    SquareMesh mesh; // Declare the mesh.
    
    // Setup the mesh from geometry:
    mesh.addRectangle(-10, 10, -5, 5, BND_MIRROR);
    mesh.addDisk(10, 5, 7.5, BND_MIRROR);
    mesh.removeDisk(10, 5, 2.5);
    
    // Setup boundary conditions (only vacuum because mirror is default):
    mesh.setBoundaryDisk(-10, 0, 2.5, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(-3, 5, 4.5, DIR_NORTH, BND_INPUT);
    mesh.setBoundaryRectangle( 10,  10, -5, 5, DIR_EAST, BND_OUTPUT);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.finalize();
    
    double kh = PI/2;
    double holscat = 0.05;
    double holabso = 0.;  // Absorption is zero.
    double density = 1.;
    
    WaveSystem sys("TestUnitarity", mesh, kh, density, holscat, holabso);
    
    return sys;
}

/**
 * Test the unitarity in the absence of losses.
 */
int testUnitarity() {
    std::cout << "====== TEST UNITARITY ======\n";
    
    WaveSystem sys = createUnitarySystem();
    sys.infoHamiltonian();
    
    // 1. Check unitarity:
    const uint64_t nseed = 3;   // Number of realizations of the disorder used to check unitarity.
    sys.checkUnitarity(false);  // First check the unitarity without disorder.
    
    for (uint64_t seed = 1; seed <= nseed; seed++) {// Loop over realizations of the disorder.
        sys.setDisorder(seed);
        sys.checkUnitarity(false);
    }
    
    // 2. Plot the mesh for checking:
    sys.plotMesh();
    
    return 0;
}

/**
 * Test the creation of reflective outputs using additional potential.
 */
int testReflectiveBarrier() {
    std::cout << "====== TEST REFLECTIVE BARRIER ======\n";
    
    SquareMesh mesh; // Declare the mesh.
    
    // Parameters of the mesh:
    const int nx = 160;
    const int ny = 100;
    const int xbar = nx;      // Position of the barrier.
    const double nprop = 1.5;   // Number of propagating modes, Nprop = kW/pi.
    const double kh = PI*nprop/ny;  // Value of k*h, deduced from Nprop: Nprop = kh*ny/pi.
    const double rbar = 0.04;   // Reflection probability of the barrier.
    const double holscat = 0.;  // Disorder strength, value of h/lscat.
    const double holabso = 0.;  // Absorption strength, value of h/labso.
    const double density = 1.;  // Density of the disorder.
    const std::string sysname = "test-barrier"; // Name of the system.
    
    // Setup the mesh from geometry:
    mesh.addRectangle(1, nx, 1, ny, BND_MIRROR);
    mesh.setBoundaryRectangle(1, 1, 1, ny, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(nx, nx, 1, ny, DIR_EAST, BND_OUTPUT);
    mesh.finalize();
    
    WaveSystem sys(sysname, mesh, kh, density, holscat, holabso);
    
    // Setup the reflective barrier:
    const int ninputprop = sys.getNInputProp();
    const double d2pkh2 = laplacianEigenvalue(0, ny) + kh*kh;  // Eigenvalue of D_y^2 + (kh)^2 for the fundamental mode.
    const double sinklh = std::sqrt( d2pkh2 * ( 1. - d2pkh2/4. ) );    // Compute sin(K_x).
    const dcomplex uh2 = 2. * std::sqrt(rbar/(1.-rbar)) * sinklh;    // Value U(x)*h^2 of the potential barrier.
    
    for (int y = 1; y <= ny; y++) {// Loop on the points to add to the potential barrier.
        sys.addPotential(xbar, y, uh2);
    }
    
    sys.infoHamiltonian();    
    sys.checkUnitarity(false);
    
    // Compute reflection eigenvalues:
    const int nrval = ninputprop;
    ComplexMatrix rmat(ninputprop, ninputprop), u(ninputprop, ninputprop), vh(ninputprop, ninputprop);
    RealMatrix rval(nrval, 1);
    sys.reflectionMatrix(rmat);
    svd(rmat, rval, u, vh);
    for (int i = 0; i < nrval; i++) {
        rval(i, 0) = rval(i, 0)*rval(i, 0);  // Convert singular values of "t" to reflection eigenvalues.
    }
    rval.transpose().print("Reflection eigenvalues");
    std::cout << TAG_INFO << "Mean reflection, Ravg = " << rval.mean() << "\n";
    sys.plotTransmissionStates(nrval);
    
    return 0;
}

/**
 * Main function of the test.
 */
int main(int argc, char** argv) {
    
    //testHamiltonian();
    //testUnitarity();
    testReflectiveBarrier();
    
    return 0;
}
