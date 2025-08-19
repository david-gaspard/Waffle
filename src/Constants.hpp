/****
 * @date Created on 2025-07-01 at 15:20:52 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing global constants.
 ***/
#ifndef _CONSTANTS_H
#define _CONSTANTS_H
#include <complex>
#include <string>

// Define the complex type:
typedef std::complex<double> dcomplex;

// Printing tags:
static const std::string TAG_INFO = "[INFO] ";                     // Information tag.
static const std::string TAG_WARN = "[\033[1;93mWARN\033[0m] ";    // Warning tag.
static const std::string TAG_ERROR = "[\033[1;31mERROR\033[0m] ";  // Error tag.
static const std::string TAG_EXEC = "[\033[1;95mEXEC\033[0m] ";    // Execution tag.
static const std::string PROGRAM_NAME_SHORT = "Waffle v0.1";       // Program name.
static const std::string PROGRAM_NAME_FULL  = PROGRAM_NAME_SHORT + " - Wave field from finite elements";
static const std::string PROGRAM_COPYRIGHT  = PROGRAM_NAME_SHORT + " (c) 2025 David GASPARD (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>";

// Mathematical constants:
static const double   PI = 3.1415926535897932384626;  // The fundamental constant Pi = 3.1415926535897932384626...
static const dcomplex I  = dcomplex(0.0, 1.0);  // Define the imaginary unit.

// Physical constants:
static const int    DIMENSION  = 2;          // Number of spatial dimensions of the mesh, which is always equal to 2 in the present program.
static const double FREEDOS    = 1./(4*PI);  // Free-space density of states in two dimensions. FREEDOS = nu = S_d k^(d-2) / 2*(2*pi)^d

// Numerical constants:
static const double MEPS    = 1.11e-16;  // Machine epsilon in double precision, MEPS = 2^(-53).
static const double SQRTEPS = 1.053e-8;  // Square root of the machine epsilon in double precision, SQRTEPS = MEPS^(1/2) = 2^(-53/2), used for discrete derivatives with forward difference formula.
static const double CBRTEPS = 4.806e-6;  // Cubic root of the machine epsilon in double precision, CBRTEPS = MPS^(1/3) = 2^(-53/3), used for discrete derivatives with central difference formula.

#endif
