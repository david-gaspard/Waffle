/****
 * @date Created on 2025-07-03 at 13:57:17 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the MeshPoint methods.
 ***/
#include "MeshPoint.hpp"

/**
 * Constructor of a MeshPoint with default unset boundaries.
 */
MeshPoint::MeshPoint(const int x, const int y, const int bndtype) {
    
    this->x = x;
    this->y = y;
    
    north = bndtype;
    south = bndtype;
    east  = bndtype;
    west  = bndtype;
}

/**
 * Default constructor of MeshPoint.
 */
MeshPoint::MeshPoint() : MeshPoint(0, 0, BND_MIRROR) {}

/**
 * Returns the neighbor index in the given direction.
 */
int MeshPoint::neighbor(const Direction dir) const {
    switch (dir) {
        case DIR_NORTH:
            return north;
        case DIR_SOUTH:
            return south;
        case DIR_EAST:
            return east;
        case DIR_WEST:
            return west;
        default:
            return -1;  // All direction cannot be opened at the same time (or the mesh only contains 1 point).
    }
}

/**
 * Returns "true" if the current point possesses an open boundary condition (either "open", "input", or "output") in the given direction.
 */
bool MeshPoint::isOpening(const Direction dir) const {
    switch (dir) {
        case DIR_NORTH:
            return north == BND_OPEN || north == BND_INPUT || north == BND_OUTPUT;
        case DIR_SOUTH:
            return south == BND_OPEN || south == BND_INPUT || south == BND_OUTPUT;
        case DIR_EAST:
            return east  == BND_OPEN || east  == BND_INPUT || east  == BND_OUTPUT;
        case DIR_WEST:
            return west  == BND_OPEN || west  == BND_INPUT || west  == BND_OUTPUT;
        default:
            return false;  // All direction cannot be opened at the same time (or the mesh only contains 1 point).
    }
}

/**
 * Returns "true" is the current point belongs to an opening, whatever the direction.
 */
bool MeshPoint::isOpening() const {
    return isOpening(DIR_NORTH) || isOpening(DIR_SOUTH) || isOpening(DIR_EAST) || isOpening(DIR_WEST);
}

/**
 * Returns the opposite direction of the given direction "dir".
 */
Direction opposite(const Direction dir) {
    switch (dir) {
        case DIR_NORTH:
            return DIR_SOUTH;
        case DIR_SOUTH:
            return DIR_NORTH;
        case DIR_EAST:
            return DIR_WEST;
        case DIR_WEST:
            return DIR_EAST;
        default:
            return DIR_ALL;
    }
}

/**
 * Returns a string version of the given direction.
 */
const std::string directionString(const Direction dir) {
    switch (dir) {
        case DIR_NORTH:
            return "north";
        case DIR_SOUTH:
            return "south";
        case DIR_EAST:
            return "east";
        case DIR_WEST:
            return "west";
        case DIR_ALL:
            return "all";
        default:
            return "direction?";
    }
}

/**
 * Returns a string version of the given boundary type.
 */
const std::string boundaryTypeString(const int bndtype) {
    switch (bndtype) {
        case BND_MIRROR:
            return "mirror";
        case BND_OPEN:
            return "open";
        case BND_INPUT:
            return "input";
        case BND_OUTPUT:
            return "output";
        default:
            return std::to_string(bndtype);
    }
}

/**
 * Comparator function used to sort the points in column-major ordering, i.e., from up to down, then left to right.
 * This function returns "true" when "p1" is located before "p2" according to the column-major order.
 */
bool comparePoints(const MeshPoint& p1, const MeshPoint& p2) {
    return p1.x < p2.x || (p1.x == p2.x && p1.y > p2.y);
}
