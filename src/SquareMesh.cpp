/****
 * @date Created on 2025-06-27 at 10:06:14 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing a two-dimensional square mesh object.
 ***/
#include "SquareMesh.hpp"
#include "BaseTools.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

/**
 * Constructor of the Square Mesh object.
 */
SquareMesh::SquareMesh() {
    //std::cout << TAG_INFO << "Creating SquareMesh...\n";
    ready = false;
}

/**
 * Destructor.
 */
SquareMesh::~SquareMesh() {
    //std::cout << TAG_INFO << "Deleting SquareMesh...\n";
}

/**
 * Check the ready status of the present SquareMesh object and stops the program if not ready.
 * This is a safety measure to avoid accessing data from a partially initialized SquareMesh object.
 */
void SquareMesh::checkReady(const std::string& name) const {
    if (not ready) {
        throw std::logic_error("In " + name + ": SquareMesh is not completely initialized. Please use finalize().");
    }
}

/**
 * Return the point at position "i" in the mesh.
 */
MeshPoint SquareMesh::getPoint(const unsigned int i) const {
    checkReady("getPoint()");
    return point.at(i);  // vector objects already make bound checking.
}

/**
 * Returns the number of points in the mesh.
 */
unsigned int SquareMesh::getNPoint() const {
    checkReady("getNPoint()");
    return point.size();
}

/**
 * Returns the number of boundary conditions of a certain type in the mesh.
 * This method can be used for instance to determine the number of input/output channels in order to compute
 * the intensity profile of transmission eigenchannels (see: saveIntensities()).
 */
unsigned int SquareMesh::getNBoundary(const int bndtype) const {
    
    checkReady("getNBoundary()");
    
    unsigned int nbnd = 0;
    
    for (MeshPoint p : point) {// Loop on the points, looking for boundary conditions.
        nbnd += (p.north == bndtype) + (p.south == bndtype) + (p.east == bndtype) + (p.west == bndtype);
    }
    
    return nbnd;
}

/**
 * Returns a copy of the "opening" vector, i.e., a list of openings, each one containing the indices of points 
 * with open boundary conditions (either "open", "input", or "output").
 */
std::vector<Opening> SquareMesh::getOpening() const {
    checkReady("getOpening()");
    return opening;
}

/**
 * Returns the number of openings in the present SquareMesh object.
 */
unsigned int SquareMesh::getNOpening() const {
    checkReady("getNOpening()");
    return opening.size();
}

/**
 * Return the index of a given point in the mesh.
 * If the point does not exist, then return BND_DEFAULT (a negative value).
 */
int SquareMesh::indexOf(const int x, const int y) const {
    MeshPoint p;
    for (unsigned int i = 0; i < point.size(); i++) {
        p = point.at(i);
        if (p.x == x && p.y == y){
            return i;
        }
    }
    return -1;
}

/**
 * Returns true if the mesh contains the given point (x, y).
 */
bool SquareMesh::containsPoint(const int x, const int y) const {
    return indexOf(x, y) >= 0;
}

/**
 * Add a single point to the mesh only if it is not already contained in the mesh.
 */
void SquareMesh::addPoint(const int x, const int y, const int bndtype) {
    if (not containsPoint(x, y)) {
        ready = false;  // If the mesh is changed, then the neighbors must be recomputed.
        point.push_back(MeshPoint(x, y, bndtype)); // Add point with default neighboring value.
    }
}

/**
 * Add a square region to the mesh.
 */
void SquareMesh::addRectangle(int xmin, int xmax, int ymin, int ymax, const int bndtype) {
    
    if (xmin > xmax) std::swap(xmin, xmax);
    if (ymin > ymax) std::swap(ymin, ymax);
    
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            addPoint(x, y, bndtype);
        }
    }
}

/**
 * Add a circular region to the mesh.
 */
void SquareMesh::addDisk(int x0, int y0, double radius, const int bndtype) {
    int xmin = (int) std::floor(x0 - radius);
    int xmax = (int) std::ceil (x0 + radius);
    int ymin = (int) std::floor(y0 - radius);
    int ymax = (int) std::ceil (y0 + radius);
    int r2 = (int) std::ceil(radius*radius);
    int dx, dy;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            dx = x - x0;
            dy = y - y0;
            if (dx*dx + dy*dy <= r2) {
                addPoint(x, y, bndtype);
            }
        }
    }
}

/**
 * Add a polygon region to the mesh. Uses the even-odd filling rule.
 */
void SquareMesh::addPolygon(const std::vector<Vector2D>& polygon, const int bndtype) {
    
    // 1. Determine the bounds of the polygon:
    int xmin, xmax, ymin, ymax;
    polygonBounds(polygon, xmin, xmax, ymin, ymax);
    
    // 2. Loop on points in the rectangular region:
    Vector2D p;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            p = Vector2D(x, y);   // Current point that we attempt to add to the mesh.
            if (p.windingNumber(polygon) % 2 != 0) {// Uses the even-odd rule to fill the polygon (fill if the winding number is odd).
                addPoint(x, y, bndtype);
            }
        }
    }
}

/**
 * Add a polygon using the coordinates given by a file.
 */
void SquareMesh::addPolygon(const char* filename, const double scale, const int bndtype) {
    
    std::vector<Vector2D> polygon;
    std::ifstream ifs(filename);  // Open the file.
    std::string line;
    double x, y;
    size_t sz;
    
    if (ifs.fail()) {// Check if the file is opened.
        std::string info = "In addPolygon(): Cannot open file '" + std::string(filename) + "'.";
        throw std::runtime_error(info);
    }
    
    while (getline(ifs, line)) {// Loop on each line of the file.
        try {
            x = scale * std::stod(line, &sz);
            y = scale * std::stod(line.substr(sz));
            polygon.push_back(Vector2D(x, y));
        }
        catch (const std::invalid_argument& exc) {
            std::cerr << TAG_WARN << "In addPolygon(): Argument is not a number, skipping line (what = '" 
                      << exc.what() << "').\n";
        }
    }
    
    ifs.close();
    
    addPolygon(polygon, bndtype);  // Add the parsed polygon.
}

/**
 * Remove the point at a given position (x, y) and assigns the boundary conditions to the neighboring points.
 */
void SquareMesh::removePoint(const int x, const int y) {
    
    int i = indexOf(x, y);  // Find the index of the point.
    
    if (i >= 0) {
        ready = false;  // If the mesh is changed, then the neighbors must be recomputed.
        point.erase(point.begin() + i);
    }
}

/**
 * Remove a rectangular region of points.
 */
void SquareMesh::removeRectangle(int xmin, int xmax, int ymin, int ymax) {
    
    if (xmin > xmax) std::swap(xmin, xmax);
    if (ymin > ymax) std::swap(ymin, ymax);
    
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            removePoint(x, y);
        }
    }
}

/**
 * Remove a circular region from the mesh.
 */
void SquareMesh::removeDisk(const int x0, const int y0, const double radius) {
    int xmin = (int) std::floor(x0 - radius);
    int xmax = (int) std::ceil (x0 + radius);
    int ymin = (int) std::floor(y0 - radius);
    int ymax = (int) std::ceil (y0 + radius);
    int r2 = (int) std::ceil(radius*radius);
    int dx, dy;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            dx = x - x0;
            dy = y - y0;
            if (dx*dx + dy*dy <= r2) {
                removePoint(x, y);
            }
        }
    }
}

/**
 * Remove half a circular region from the mesh (only the upper half disk).
 */
void SquareMesh::removeHalfDisk(const int x0, const int y0, const double radius) {
    int xmin = (int) std::floor(x0 - radius);
    int xmax = (int) std::ceil (x0 + radius);
    int ymin = y0;
    int ymax = (int) std::ceil (y0 + radius);
    int r2 = (int) std::ceil(radius*radius);
    int dx, dy;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            dx = x - x0;
            dy = y - y0;
            if (dx*dx + dy*dy <= r2) {
                removePoint(x, y);
            }
        }
    }
}

/**
 * Removes a polygon region from the mesh. Uses the even-odd rule.
 */
void SquareMesh::removePolygon(const std::vector<Vector2D>& polygon) {
    
    // 1. Determine the bounds of the polygon:
    int xmin, xmax, ymin, ymax;
    polygonBounds(polygon, xmin, xmax, ymin, ymax);
    
    // 2. Loop on points in the rectangular region:
    Vector2D p;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            p = Vector2D(x, y);   // Current point that we attempt to add to the mesh.
            if (p.windingNumber(polygon) % 2 != 0) {// Uses the even-odd rule to remove the polygon.
                removePoint(x, y);
            }
        }
    }
}

/**
 * Setup the boundary condition for the point at the given position (x, y).
 */
void SquareMesh::setBoundaryPoint(const int x, const int y, const Direction dir, const int bndtype) {
    int i = indexOf(x, y);
    if (i >= 0) {// Only works if the point is in the mesh.
        ready = false;  // If the boundary conditions are changed, then the openings must be recomputed.
        MeshPoint& p = point.at(i);
        if ( (dir == DIR_ALL || dir == DIR_NORTH) && p.north < 0) p.north = bndtype;
        if ( (dir == DIR_ALL || dir == DIR_SOUTH) && p.south < 0) p.south = bndtype;
        if ( (dir == DIR_ALL || dir == DIR_EAST ) && p.east  < 0) p.east  = bndtype;
        if ( (dir == DIR_ALL || dir == DIR_WEST ) && p.west  < 0) p.west  = bndtype;
    }
}

/**
 * Setup the same boundary condition for a rectangular region of points.
 */
void SquareMesh::setBoundaryRectangle(int xmin, int xmax, int ymin, int ymax, const Direction dir, const int bndtype) {
    
    if (xmin > xmax) std::swap(xmin, xmax);
    if (ymin > ymax) std::swap(ymin, ymax);
    
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            setBoundaryPoint(x, y, dir, bndtype);
        }
    }
}

/**
 * Setup the boundary condition 'bndtype' in a disk region of points.
 */
void SquareMesh::setBoundaryDisk(const int x0, const int y0, const double radius, const Direction dir, const int bndtype) {
    int xmin = (int) std::floor(x0 - radius);
    int xmax = (int)  std::ceil(x0 + radius);
    int ymin = (int) std::floor(y0 - radius);
    int ymax = (int)  std::ceil(y0 + radius);
    int r2 = (int) std::ceil(radius*radius);
    int dx, dy;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            dx = x - x0;
            dy = y - y0;
            if (dx*dx + dy*dy <= r2) {
                setBoundaryPoint(x, y, dir, bndtype);
            }
        }
    }
}

/**
 * Sort the points of the mesh in reading order, i.e., from up to down, from left to right.
 * The purpose of this method is to make the identification of ports easier [see also: computePorts()].
 * This function is also possibly helpful to place the Hamiltonian matrix elements as close to the diagonal as possible, and
 * is therefore advantageous when inverting the Hamiltonian matrix.
 */
void SquareMesh::sortPoints() {
    std::sort(point.begin(), point.end(), comparePoints);
}

/**
 * Identify the nearest neighbors of each point in the mesh.
 */
void SquareMesh::fixNeighbors() {
    int inorth, isouth, ieast, iwest;
    for (MeshPoint& p : point) {// Loop over the points. Note the address because the point is modified.
        inorth = indexOf(p.x, p.y+1);
        isouth = indexOf(p.x, p.y-1);
        ieast  = indexOf(p.x+1, p.y);
        iwest  = indexOf(p.x-1, p.y);
        if (inorth >= 0) p.north = inorth;
        if (isouth >= 0) p.south = isouth;
        if (ieast  >= 0) p.east  = ieast;
        if (iwest  >= 0) p.west  = iwest;
    }
}

/**
 * Prepare the "ports" by identifying and classifying the points with open, input, and output boundary conditions.
 * Note that a "port" is, by definition, an alignment of nearest-neighboring points with the same boundary condition
 * (either "open", "input", or "output") in the same direction.
 * Warning: This method requires the points to be sorted in reading order to operate properly because it is essential
 * that neighboring points are traveled successively in a determined order.
 */
void SquareMesh::computeOpening() {
    
    MeshPoint p, lastp;
    bool added;
    
    opening.clear();  // Clears the opening in case it was already computed before.
    
    for (unsigned int i = 0; i < point.size(); i++) {// Loop over the points using their index.
        p = point.at(i); // Extract the point
        for (Direction dir : allDirections) {// Loop over all the directions.
            if (p.isOpening(dir)) {// If the point contains an open boundary condition ("open", "input", or "output").
                added = false; // By default the point is not added to the list of openings.
                for (Opening& op : opening) {// Loop over the list of openings. Note the address because "op" is modified.
                    
                    lastp = point.at(op.index.back());  // Extract the last point of opening "iop".
                    
                    if (op.direction == dir && op.bndtype == p.neighbor(dir) && std::abs(p.x - lastp.x) + std::abs(p.y - lastp.y) == 1) {
                            // If "op" is an opening in the same direction as "p", if "op" has the same boundary condition as "p", 
                            // and if "lastp" is the nearest neighbor of the current point "p".
                        op.index.push_back(i);  // Add the index to the opening "iop".
                        added = true;  // Signal the point as being added.
                    }
                }
                if (not added) {// If "p" not added to the list of existing openings, then create a new opening.
                    Opening newop = {};
                    newop.index.push_back(i);
                    newop.direction = dir;
                    newop.bndtype = p.neighbor(dir);
                    opening.push_back(newop);
                }
            }
        }
    }
}

/**
 * Finalize the mesh by sorting the points, fixing the neighbors, and computing the ports.
 * This method must be called after the mesh construction methods add*(), remove*(), and setBoundary*().
 */
void SquareMesh::finalize() {
    //std::cout << TAG_INFO << "Finalizing the mesh...\n";
    sortPoints();
    fixNeighbors();
    computeOpening();
    ready = true; // The SquareMesh gets ready for computations.
}

/**
 * Print a summary of the mesh.
 */
void SquareMesh::printSummary() const {
    std::cout << TAG_INFO << "SquareMesh with Npoint=" << point.size() << ", Nopening=" << opening.size() << ".\n";
}

/**
 * Print the openings to standard output. This function is mainly used for testing purposes.
 */
void SquareMesh::printOpening() const {
    
    checkReady("printOpening()");
    
    MeshPoint p;
    
    for (unsigned int iop = 0; iop < opening.size(); iop++) {// Loop over openings.
        
        std::cout << TAG_INFO << "Opening #" << iop << ", bnd=" << boundaryTypeString(opening.at(iop).bndtype)
                  << ", dir=" << directionString(opening.at(iop).direction) << " : [\n";
        
        for (unsigned int idx : opening.at(iop).index) {// Loop over the point indices composing the opening.
            p = point.at(idx); // Corresponding point in the mesh.
            std::cout << "\t" << idx << " (" << p.x << ", " << p.y << ")" 
                      << ", north=" << boundaryTypeString(p.north) << ", south=" << boundaryTypeString(p.south) 
                      << ", east=" << boundaryTypeString(p.east) << ", west=" << boundaryTypeString(p.west) << "\n";
        }
        
        std::cout << "]\n";
    }
    
}

/**
 * Save the points to a file "filename" using the separator "sep".
 */
void SquareMesh::saveMesh(const std::string& filename, const char* sep) const {
    
    checkReady("saveMesh()");
    
    std::ofstream ofs;
    ofs.open(filename.c_str());
    
    writeTimestamp(ofs, "%% ");
    
    ofs << "%% SquareMesh with Npoint=" << point.size() << ".\n"
        << "x" << sep << "y" << sep << "north" << sep << "south" << sep << "east" << sep << "west\n";
    
    for (MeshPoint p : point) {
        ofs << p.x << sep << p.y << sep << boundaryTypeString(p.north) << sep 
                                        << boundaryTypeString(p.south) << sep 
                                        << boundaryTypeString(p.east) << sep 
                                        << boundaryTypeString(p.west) << "\n";
    }
    
    ofs.close();
}

/**
 * Save the points to a file "filename" using the separator "sep".
 * This is a short version for quick plots with less post-processing.
 */
void SquareMesh::saveMeshShort(const std::string& filename, const char* sep) const {
    
    checkReady("saveMeshShort()");
    
    std::ofstream ofs;
    ofs.open(filename.c_str());
    
    writeTimestamp(ofs, "%% ");
    
    ofs << "%% SquareMesh with Npoint=" << point.size() << ".\n" 
        << "x" << sep << "y" << sep << "bnd\n";
    
    for (MeshPoint p : point) {
        ofs << p.x << sep << p.y << sep;
        if (p.north == BND_INPUT || p.south == BND_INPUT || p.east == BND_INPUT || p.west == BND_INPUT) {
            ofs << "input";
        }
        else if (p.north == BND_OUTPUT || p.south == BND_OUTPUT || p.east == BND_OUTPUT || p.west == BND_OUTPUT) {
            ofs << "output";
        }
        else if (p.north == BND_OPEN || p.south == BND_OPEN || p.east == BND_OPEN || p.west == BND_OPEN) {
            ofs << "open";
        }
        else if (p.north == BND_MIRROR || p.south == BND_MIRROR || p.east == BND_MIRROR || p.west == BND_MIRROR) {
            ofs << "mirror";
        }
        else {
            ofs << "bulk";
        }
        ofs << "\n";
    }
    
    ofs.close();
}

/**
 * Plot the mesh using external plotting script.
 * The name of the CSV file written by saveMesh() is "filename".
 */
void SquareMesh::plotMesh(const std::string& filename) const {
    saveMesh(filename, ", ");
    std::string cmd("plot/plot_mesh.py " + filename);
    std::cout << TAG_EXEC << cmd << "\n";
    if (std::system(cmd.c_str())) {
        std::cout << TAG_WARN << "The plot script returned an error.\n";
    }
}
