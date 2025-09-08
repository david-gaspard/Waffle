/****
 * @date Created on 2025-09-05 at 19:15:57 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing color conversion tools.
 ***/
#include "Color.hpp"

/**************************************************************************************************
 * COLOR MANAGEMENT AND CONVERSION
 **************************************************************************************************/

/**
 * Convert the given RGB primaries in the interval [0, 1] to standard 8-bit RGB color.
 */
Color RGBColor(double red, double green, double blue) {
    
    // Clip the values to the interval [0, 1] to avoid exceeding the bounds:
    if      (red   < 0.) red   = 0.;
    else if (red   > 1.) red   = 1.;
    if      (green < 0.) green = 0.;
    else if (green > 1.) green = 1.;
    if      (blue  < 0.) blue  = 0.;
    else if (blue  > 1.) blue  = 1.;
    
    return Color{(uint8_t)std::round(MAX_COLOR*red), (uint8_t)std::round(MAX_COLOR*green), (uint8_t)std::round(MAX_COLOR*blue)};
}

/**
 * Convert the given HSV (hue/saturation/value) values defined on the interval [0, 1] to standard 8-bit RGB color.
 * Conventions are such that hue=0 and hue=1 are the red color, and hue=0.5 is the cyan color.
 * Saturation "sat" controls the amount of white, and "val" the amount of black.
 */
Color HSVColor(double hue, double sat, double val) {
    
    // Clip the values to avoid exceeding the bounds:
    if      (hue < 0.) hue = 0.;
    else if (hue > 1.) hue = 1.;
    if      (sat < 0.) sat = 0.;
    else if (sat > 1.) sat = 1.;
    if      (val < 0.) val = 0.;
    else if (val > 1.) val = 1.;
    
    // Convert the hue to RGB values:
    hue = 6.*hue;                         // hue is now in the interval [0, 6].
    int hue_int = (int)std::floor(hue);   // hue_int is an integer {0, 1, 2, 3, 4, 5}.
    double hue_frac = hue - hue_int;      // hue_frac is a double in the interval [0, 1].
    double r, g, b;  // RGB values in the interval [0, 1].
    
    switch (hue_int) {
        case 0:
            r = 1; g = hue_frac; b = 0;
            break;
        case 1:
            r = 1 - hue_frac; g = 1; b = 0;
            break;
        case 2:
            r = 0; g = 1; b = hue_frac;
            break;
        case 3:
            r = 0; g = 1 - hue_frac; b = 1;
            break;
        case 4:
            r = hue_frac; g = 0; b = 1;
            break;
        case 5:
            r = 1; g = 0; b = 1 - hue_frac;
            break;
        default:
            r = 1; g = hue_frac; b = 0;
    }
    
    return RGBColor(val*((1. - sat) + sat*r), val*((1. - sat) + sat*g), val*((1. - sat) + sat*b));
}

/**
 * Returns a color representation of the given complex number.
 * First version using only the hue (full saturation, full value).
 */
Color complexColor1(const dcomplex z) {
    
    double hue = std::arg(z)/(2*PI);  // Hue is in the interval [-0.5, 0.5].
    if (hue < 0.) hue += 1.;  // Hue is now in [0, 1].
    
    return HSVColor(hue, 1., 1.);
}

/**
 * Returns a color representation of the given complex number.
 * Second version using hue but also fading to white for small values and to black for large values.
 * Full saturation and full value occur only for unit complex numbers (|z|=1).
 */
Color complexColor2(const dcomplex z) {
    
    double hue = std::arg(z)/(2*PI);  // Hue is in the interval [-0.5, 0.5].
    if (hue < 0.) hue += 1.;  // Hue is now in [0, 1].
    
    double sat = (4./PI) * std::atan(std::abs(z));
    double val = 2. - sat;
    
    return HSVColor(hue, sat, val);
}

/**
 * Returns true if the two colors are equal, false otherwise.
 */
bool operator==(const Color& c1, const Color& c2) {
    return c1.r == c2.r && c1.g == c2.g && c1.b == c2.b;
}

/**
 * Overload the stream operator to print the color in a portable pixmap file, a PPM file (see: https://en.wikipedia.org/wiki/Netpbm).
 */
std::ofstream& operator<<(std::ofstream& ofs, const Color& c) {
    ofs << c.r << c.g << c.b;
    return ofs;
}
