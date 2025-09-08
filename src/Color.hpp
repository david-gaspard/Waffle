/****
 * @date Created on 2025-09-05 at 19:12:40 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the color structure and color conversion tools.
 ***/
#ifndef _COLOR_H
#define _COLOR_H
#include "Constants.hpp"
#include <fstream>

/**
 * Structure for a RGB color (8-bit depth).
 */
struct Color { uint8_t r; uint8_t g; uint8_t b; };

static const uint8_t MAX_COLOR = 255;  // Maximum color value (8-bit depth).

/**
 * Standard colors:
 */
static const Color COLOR_BLACK{0, 0, 0};
static const Color COLOR_WHITE{MAX_COLOR, MAX_COLOR, MAX_COLOR};
static const Color COLOR_RED{MAX_COLOR, 0, 0};
static const Color COLOR_GREEN{0, MAX_COLOR, 0};
static const Color COLOR_BLUE{0, 0, MAX_COLOR};
static const Color COLOR_YELLOW{MAX_COLOR, MAX_COLOR, 0};
static const Color COLOR_CYAN{0, MAX_COLOR, MAX_COLOR};
static const Color COLOR_MAGENTA{MAX_COLOR, 0, MAX_COLOR};

/**
 * Conversion functions:
 */
Color RGBColor(double red, double green, double blue);
Color HSVColor(double hue, double sat, double val);
Color complexColor1(const dcomplex z);
Color complexColor2(const dcomplex z);

bool operator==(const Color& c1, const Color& c2);
std::ofstream& operator<<(std::ofstream& ofs, const Color& c);

#endif
