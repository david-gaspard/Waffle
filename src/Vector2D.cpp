/****
 * @date Created on 2025-07-03 at 13:38:38 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the Vector2D methods.
 ***/
#include "Vector2D.hpp"
#include <cmath>

/**
 * Constructor of the Vector2D object.
 */
Vector2D::Vector2D(int x, int y) : x((double)x), y((double)y) {}

/**
 * Constructor of the Vector2D structure.
 */
Vector2D::Vector2D(double x, double y) : x(x), y(y) {}

/**
 * Constructor of the Vector2D structure from a mesh point.
 */
Vector2D::Vector2D(MeshPoint p) : x((double)p.x), y((double)p.y) {}

/**
 * Default constructor of the Vector2D structure.
 */
Vector2D::Vector2D() : x(0.), y(0.) {}

/**
 * Addition of two vectors.
 */
Vector2D Vector2D::operator+(const Vector2D& v) const {
    return Vector2D(x + v.x, y + v.y);
}

/**
 * Subtraction of two vectors.
 */
Vector2D Vector2D::operator-(const Vector2D& v) const {
    return Vector2D(x - v.x, y - v.y);
}

/**
 * Norm of a vector.
 */
double Vector2D::norm() const {
    return sqrt(x*x + y*y);
}

/**
 * Return the normalized vector.
 */
Vector2D Vector2D::normalize() const {
    double n = norm();
    return Vector2D(x/n, y/n);
}

/**
 * Dot (scalar) product of two vectors.
 */
double Vector2D::dot(const Vector2D& v) const {
    return x * v.x + y * v.y;
}

/**
 * Cross (pseudo-scalar) product of two vectors.
 */
double Vector2D::cross(const Vector2D& v) const {
    return x * v.y - y * v.x;
}

/**
 * Compute the bounds of a polygon. std::vector will return an out of bound error if the polygon vector is empty.
 */
void polygonBounds(const std::vector<Vector2D>& polygon, int& xmin, int& xmax, int& ymin, int& ymax) {
    
    double px, py;
    int fpx, cpx, fpy, cpy;
    
    px = polygon.at(0).x;  // Extract the (x, y) coordinates of the first vertex.
    py = polygon.at(0).y;  // Note: std::vector will return an out of bound error if the polygon vector is empty.
    
    xmin = (int) std::floor(px);  // Initialize the minima and maxima.
    xmax = (int)  std::ceil(px);
    ymin = (int) std::floor(py);
    ymax = (int)  std::ceil(py);
    
    for (unsigned int i = 1; i < polygon.size(); i++) {// Loop on the vertices of the polygon.
        
        px = polygon.at(i).x;  // Extract the (x, y) coordinates of the current vertex.
        py = polygon.at(i).y;
        
        fpx = (int) std::floor(px);
        cpx = (int)  std::ceil(px);
        fpy = (int) std::floor(py);
        cpy = (int)  std::ceil(py);
        
        if (fpx < xmin) xmin = fpx;
        if (cpx > xmax) xmax = cpx;
        if (fpy < ymin) ymin = fpy;
        if (cpy > ymax) ymax = cpy;
    }
}

/**
 * Returns the oriented winding number of the given polygon with respect to the present point Vector2D.
 * The winding number is the number of turns the polygon does around a point. It is positive for counterclockwise direction, and negative otherwise.
 * 
 * The winding number can be used to determine if a point is inside a polygon depending on the filling rule:
 * - The nonzero rule: inside = (windingNumber() != 0)
 * - The evenodd rule: inside = (windingNumber()%2 != 0)
 * 
 * The algorithm is inspired by: K. Hormann and A. Agathos, Comput. Geom. 20, 131-144 (2001), https://doi.org/10.1016/S0925-7721(01)00012-8
 * See also: https://en.wikipedia.org/wiki/Point_in_polygon
 * See also: https://en.wikipedia.org/wiki/Even-odd_rule
 */
int Vector2D::windingNumber(const std::vector<Vector2D>& polygon) const {
    
    int w = 0;  // Initialize the winding number to zero at the beginning.
    
    double p0x, p0y, p1x, p1y, c;
    int i0, n = polygon.size();  // Number of vertices (or edges) of the polygon.
    
    for (int i = 0; i < n; i++) {// Loop on the edges of the polygon. We consider the edge: (i-1) -> i.
        
        i0 = (i == 0) ? n-1 : i-1;  // If i=0, then the previous point is the last one (because the polygon is cyclic).
        p0x = polygon.at(i0).x;
        p0y = polygon.at(i0).y;
        p1x = polygon.at(i).x;
        p1y = polygon.at(i).y;
        
        if ( (p0y < y) != (p1y < y) ) {// If the path P0 -> P1 crosses the horizontal line y.
            
            c = (p0x - x) * (p1y - y) - (p0y - y) * (p1x - x); // Cross product (P0 - A) x (P1 - A). If c>0, then we rotate counterclockwise.
            
            if (c > 0.) w++; // If the crossing is counterclockwise oriented, then increases the winding number.
            else if (c < 0.) w--; // Otherwise, then decreases the winding number. Do not change for perfect intersection.
            
        }
    }
    
    return w/2;  // Return half the total winding because there are two possible crossings of the horizontal line: To the left and to the right.
}

/**
 * Overloads the stream operator for printing the vector.
 */
std::ostream& operator<<(std::ostream& os, const Vector2D& v) {
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}
