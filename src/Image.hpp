/****
 * @date Created on 2025-08-12 at 11:05:45 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing raster image utilities.
 ***/
#ifndef _IMAGE_H
#define _IMAGE_H
#include "Color.hpp"

/**
 * Class defining a raster image.
 */
class Image {
    
    private:
    
    int height;   // Height of the image. Number of pixels in vertical direction.
    int width;    // Width of the image. Number of pixels in horizontal direction.
    Color* data;  // Array of pixels in row-major format (standard reading order).
    
    public:
    
    // Constructors/Destructors:
    Image(const int height, const int width);
    Image(const std::string& filename);
    ~Image();
    
    // Getters/Setters:
    int getHeight() const;
    int getWidth() const;
    Color& operator()(const int i, const int j) const;
    
    // Export:
    void savePPM(const std::string& filename) const;
    void savePNG(const std::string& filename) const;
    
    // Image transformations:
    void fill(const Color color);
    void invert();
    
};

#endif
