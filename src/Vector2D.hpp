/****
 * @date Created on 2025-07-03 at 13:33:27 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the Vector2D structure.
 ***/
#ifndef _VECTOR2D_H
#define _VECTOR2D_H
#include "MeshPoint.hpp"
#include <vector>
#include <iostream>

/**
 * Defines a 2D vector in double precision.
 */
struct Vector2D {
    
    double x; // X coordinate of the 2D vector.
    double y; // X coordinate of the 2D vector.
    
    // Constructors:
    Vector2D(int x, int y);
    Vector2D(double x, double y);
    Vector2D(MeshPoint p);
    Vector2D();
    
    // Methods:
    Vector2D operator+(const Vector2D& v) const;
    Vector2D operator-(const Vector2D& v) const;
    double norm() const;
    Vector2D normalize() const;
    double dot(const Vector2D& v) const;
    double cross(const Vector2D& v) const;
    int windingNumber(const std::vector<Vector2D>& polygon) const;
    
    friend void polygonBounds(const std::vector<Vector2D>& polygon, int& xmin, int& xmax, int& ymin, int& ymax);
    friend std::ostream& operator<<(std::ostream& os, const Vector2D& v);
    
};

#endif
