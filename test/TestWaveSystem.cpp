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
    
    WaveSystem sys("TestHamiltonian", mesh, kh, holscat, holabso);
    
    return sys;
}

/**
 * Test the construction of the Hamiltonian.
 */
int testHamiltonian() {
    std::cout << "====== TEST HAMILTONIAN AND INDEPENDENT TERMS ======\n";
    
    WaveSystem sys = createTestSystem();
    
    sys.infoHamiltonian();  // Summary of the Hamiltonian to stdout.
    sys.plotHamiltonian();  // Plot the Hamiltonian matrix to an image file.
    sys.plotInputState();   // Plot the input state matrix to an image file.
    sys.plotOutputState();  // Plot the output state matrix to a file.
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
    
    WaveSystem sys("TestUnitarity", mesh, kh, holscat, holabso);
    
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
 * Main function of the test.
 */
int main(int argc, char** argv) {
    
    //testHamiltonian();
    testUnitarity();
    
    return 0;
}
