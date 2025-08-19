/****
 * @date Created on 2025-07-03 at 13:52:26 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the MeshPoint definition.
 ***/
#ifndef _MESHPOINT_H
#define _MESHPOINT_H
#include <string>

/**
 * Enumeration specifying the four possible directions on a square lattice (plus the omni-directional element):
 */
enum Direction {DIR_NORTH, DIR_SOUTH, DIR_EAST, DIR_WEST, DIR_ALL};

static const Direction allDirections[] = {DIR_NORTH, DIR_SOUTH, DIR_EAST, DIR_WEST};  // List of all possible directions.

/**
 * Different types of boundary conditions:
 */
static const int BND_MIRROR  = -1;  // Index of neighboring points used for mirror boundary conditions (zero current condition).
static const int BND_OPEN    = -2;  // Index of neighboring points used for open boundary conditions (extrapolated-length boundary conditions).
static const int BND_INPUT   = -3;  // Index of neighboring points used for input boundary condition.
static const int BND_OUTPUT  = -4;  // Index of neighboring points used for input boundary condition.

/**
 * Define the integer point type used in the square mesh.
 */
struct MeshPoint {
    
    int x;     // X coordinate of the point.
    int y;     // Y coordinate of the point.
    int north; // Index of the nearest neighboring point at the north (x, y+1). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    int south; // Index of the nearest neighboring point at the south (x, y-1). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    int east;  // Index of the nearest neighboring point at the east (x+1, y). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    int west;  // Index of the nearest neighboring point at the west (x-1, y). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    
    // Constructors:
    MeshPoint(const int x, const int y, const int bndtype);
    MeshPoint();
    
    // Methods:
    int neighbor(const Direction dir) const;
    bool isOpening(const Direction dir) const;
    bool isOpening() const;
    
};

Direction opposite(const Direction dir);
const std::string directionString(const Direction dir);
const std::string boundaryTypeString(const int bndtype);
bool comparePoints(const MeshPoint& p1, const MeshPoint& p2);

#endif
