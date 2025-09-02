/****
 * @date Created on 2025-06-27 at 16:23:02 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the SquareMesh class.
 ***/
#ifndef _SQUARE_MESH_H
#define _SQUARE_MESH_H
#include "Constants.hpp"
#include "MeshPoint.hpp"
#include "Vector2D.hpp"
#include <chrono>

/**
 * Defines an opening, i.e., an alignment of nearest-neighboring points with an open boundary condition ("open", "input", or "output") in the same direction.
 */
struct Opening {
    
    std::vector<unsigned int> index;  // Indices of points involved in the present opening. The indices correspond to the vector "point" (see: SquareMesh).
                                      // Note: As an alternative, this object could also store a copy of the mesh points, instead of the indices,
                                      // but then the mesh points must store their own index.
    Direction direction;  // Boundary direction of the present opening.
    int bndtype;          // Type of boundary condition (either "input", "output", or "open").
    
};

/**
 * Class defining the Square Mesh object.
 */
class SquareMesh {
    
    private:
    
    std::vector<MeshPoint> point;  // List of positions of points (belonging to the integer lattice). They also contain neighbor indices.
    std::vector<Opening> opening;  // List of openings. Each opening contains a list of indices, and a boundary direction. 
                                   // No distinction between "open", "input", and "output" boundary conditions, i.e., they are gathered in the same opening.
    
    bool ready;  // True when the SquareMesh is ready for computations. False, otherwise. This boolean becomes True after finalize() is called.
    std::chrono::time_point<std::chrono::steady_clock> start_build; // Time point at which the square mesh is built.
    
    public:
    
    // Constructors/Destructors:
    SquareMesh();
    ~SquareMesh();
    
    // Getters:
    MeshPoint getPoint(const unsigned int i) const;
    unsigned int getNPoint() const;
    unsigned int getNBoundary(const int bndtype) const;
    std::vector<Opening> getOpening() const;
    unsigned int getNOpening() const;
    
    // Determine the index of a given point:
    int indexOf(const int x, const int y) const;
    bool containsPoint(const int x, const int y) const;
    
    // Add points:
    void addPoint(const int x, const int y, const int bndtype);
    void addRectangle(int xmin, int xmax, int ymin, int ymax, const int bndtype);
    void addDisk(const int x0, const int y0, const double radius, const int bndtype);
    void addPolygon(const std::vector<Vector2D>& polygon, const int bndtype);
    void addPolygon(const char* filename, const double scale, const int bndtype);
    
    // Remove points:
    void removePoint(const int x, const int y);
    void removeRectangle(int xmin, int xmax, int ymin, int ymax);
    void removeDisk(const int x0, const int y0, const double radius);
    void removeHalfDisk(const int x0, const int y0, const double radius);
    void removePolygon(const std::vector<Vector2D>& polygon);
    
    // Assign boundary conditions:
    void setBoundaryPoint(const int x, const int y, const Direction dir, const int bndtype);
    void setBoundaryRectangle(int xmin, int xmax, int ymin, int ymax, const Direction dir, const int bndtype);
    void setBoundaryDisk(const int x0, const int y0, const double radius, const Direction dir, const int bndtype);
    
    // Finalize the mesh:
    void finalize();
    
    // Output methods:
    void printSummary() const;
    void printOpening() const;
    void saveMesh(const std::string& filename, const char* sep) const;
    void saveMeshShort(const std::string& filename, const char* sep) const;
    void plotMesh(const std::string& filename) const;
    
    // Private computational methods:
    private:
    
    void checkReady(const std::string& name) const;
    void sortPoints();
    void fixNeighbors();
    void computeOpening();
    
};

#endif
