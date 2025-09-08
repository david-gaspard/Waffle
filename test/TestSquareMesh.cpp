/****
 * @date Created on 2025-06-27 at 16:32:13 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code to test the SquareMesh class.
 ***/
#include "SquareMesh.hpp"
#include <iostream>

/**
 * Test the basic functions of SquareMesh object.
 */
int testMeshBasic() {
    std::cout << "====== TEST BASIC SQUARE MESH ======\n";
    
    SquareMesh mesh; // Declare the mesh.
    
    // Setup the mesh from geometry:
    mesh.addRectangle(-10, 10, -5, 5, BND_MIRROR);
    mesh.addDisk(10, 5, 7.5, BND_MIRROR);
    mesh.removeDisk(10, 5, 2.5);
    
    // Setup boundary conditions (only vacuum because mirror is default):
    mesh.setBoundaryRectangle(-10, -10, -5, 5, DIR_WEST, BND_OPEN);
    mesh.setBoundaryDisk(-10, 0, 2.5, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 10,  10, -5, 5, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryDisk(10, 5, 3.5, DIR_ALL, BND_OPEN);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.finalize();
    
    mesh.printOpening(); // Display the openings in terminal.
    
    // Save the mesh to a file for checking:
    //const std::string filename_short = "old/test/mesh/test_mesh_short.csv";
    //std::cout << TAG_INFO << "Saving SquareMesh to files '" << filename_short << "' and '" << filename_full << "'...\n";
    //mesh.saveMeshShort(filename_short, ", ");
    
    const std::string filename = "old/test/mesh/test_mesh.csv";
    mesh.plotMesh(filename);
    
    return 0;
}

/**
 * Test the polygon generation of SquareMesh object.
 */
int testMeshPolygon1() {
    std::cout << "====== TEST SQUARE MESH POLYGON #1 ======\n";
    
    std::vector<Vector2D> polygon = {Vector2D(10, 0), Vector2D(4, 8), Vector2D(-7, 2), Vector2D(-5, -11), Vector2D(1, -3)};
    
    SquareMesh mesh; // Declare the mesh.
    
    mesh.addPolygon(polygon, BND_MIRROR);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.finalize();
    
    // Save the mesh to a file for checking:
    const std::string filename = "old/test/mesh/test_mesh_polygon_1.csv";
    std::cout << TAG_INFO << "Saving SquareMesh to file '" << filename << "'...\n";
    mesh.saveMeshShort(filename, ", ");
    
    return 0;
}

/**
 * Test the polygon generation of SquareMesh object.
 */
int testMeshPolygon2() {
    std::cout << "====== TEST SQUARE MESH POLYGON #2 ======\n";
    
    std::vector<Vector2D> polygon = {
        Vector2D( 25, 0  ), Vector2D( 20, 12 ), Vector2D( 10, 20 ), Vector2D( -5, 18 ), Vector2D(-16, 13 ),
        Vector2D(-27, 4  ), Vector2D(-25, -6 ), Vector2D(-19, -8 ), Vector2D(-12, -15), Vector2D( -3, -32),
        Vector2D(  5, -35), Vector2D( 14, -23), Vector2D( 10, -11), Vector2D(  2, -9 ), Vector2D( -3, -15),
        Vector2D(  1, -24), Vector2D( 12, -21), Vector2D( 22, -10)
    };
    
    SquareMesh mesh; // Declare the mesh.
    
    mesh.addPolygon(polygon, BND_MIRROR);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.finalize();
    
    // Save the mesh to a file for checking:
    const std::string filename = "out/test/mesh/test_mesh_polygon_2.csv";
    std::cout << TAG_INFO << "Saving SquareMesh to file '" << filename << "'...\n";
    mesh.saveMeshShort(filename, ", ");
    
    return 0;
}

/**
 * Test the polygon generation of SquareMesh object, but using the coordinates from a file.
 */
int testMeshPolygon3() {
    std::cout << "====== TEST SQUARE MESH POLYGON #3 ======\n";
    
    SquareMesh mesh; // Declare the mesh.
    
    double scale = 60; // Scales the polygon.
    
    mesh.addPolygon("shape/eiffel-tower.csv", scale, BND_MIRROR);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.finalize();
    
    // Save the mesh to a file for checking:
    const std::string filename = "out/test/mesh/test_mesh_polygon_3.csv";
    std::cout << TAG_INFO << "Saving SquareMesh to file '" << filename << "'...\n";
    mesh.saveMeshShort(filename, ", ");
    
    return 0;
}

/**
 * Test the creation of a mesh from a PNG image.
 */
int testMeshImage() {
    
    const std::string infile = "old/test/mesh_png/mesh_1.png";
    SquareMesh mesh(infile);
    
    const std::string outfile = "old/test/mesh_png/mesh_1.csv";
    mesh.plotMesh(outfile);
    
    return 0;
}

/**
 * Main function of the test of the SquareMesh object.
 */
int main(int argc, char** argv) {
    
    //testMeshBasic();
    //testMeshPolygon1();
    //testMeshPolygon2();
    //testMeshPolygon3();
    testMeshImage();
}
