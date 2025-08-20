/****
 * @date Created on 2025-08-20 at 14:46:57 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the complete elliptic integral K(m). Note that "m" must be negative in the program and, unfortunately, this is not
 * available from the standard C++ header "cmath".
 ***/
#include "BaseTools.hpp"
#include "Constants.hpp"
#include <cmath>
#include <vector>
#include <iostream>

/**
 * Test the complete elliptic integral from the standard C++ "cmath" library.
 * Note that the standard only implements K(k) = K(sqrt(m)), i.e., the "bad" convention, so that negative values of m=k^2 cannot be reached.
 * Taken from: https://en.cppreference.com/w/cpp/numeric/special_functions/comp_ellint_1.html
 */
void testStdEllipticK() {
    std::cout << "====== TEST STD ELLIPTIC K ======\n";
    
    constexpr double π{3.1415926535897932384626};
    
    std::cout << "K(0) ≈ " << std::comp_ellint_1(0) << '\n'
              << "π/2 ≈ " << π / 2 << '\n'
              << "K(0.5) ≈ " << std::comp_ellint_1(0.5) << '\n'
              << "F(0.5, π/2) ≈ " << std::ellint_1(0.5, π / 2) << '\n'
              << "The period of a pendulum length 1m at 10° initial angle ≈ "
              << 4 * std::sqrt(1 / 9.80665) * std::comp_ellint_1(std::sin(π / 18 / 2))
              << "s,\n" "whereas the linear approximation gives ≈ "
              << 2 * π * std::sqrt(1 / 9.80665) << '\n';
}

/**
 * Test our own implementation of the complete elliptic integral K(m) using the standard convention from Wolfram Mathematica
 * which allows for negative arguments.
 */
void testOwnEllipticK() {
    std::cout << "====== TEST OWN ELLIPTIC K ======\n";
    
    // Create a checkup table of couples [m, K(m)], computed with Wolfram Mathematica :
    std::vector<std::pair<double, double>> kval = {
        {   0.99, 3.69563736298987468},
        {   0.90, 2.57809211334817319},
        {   0.50, 1.85407467730137192},
        {   0.00, 1.57079632679489662},
        {  -0.50, 1.41573720842595620},
        {  -1.00, 1.31102877714605991},
        {  -2.00, 1.17142008414676986},
        {  -5.00, 0.955503927064043934},
        { -10.00, 0.790871890238738475},
        { -20.00, 0.639782175897459880},
        { -50.00, 0.471034245408733317},
        {-100.00, 0.368219248609141033}
    };
    
    double m, km, km_expc;
    
    for (std::pair<double, double> pair : kval) {// Loop over the pairs.
        m = pair.first;
        km_expc = pair.second;
        km = ellipticK(m);
        std::cout << TAG_INFO << "K(" << m << ") \t= " << km << ",\terror = " << std::abs((km - km_expc)/km_expc) << "\n";
    }
    
}

/**
 * Main function of the test.
 */
int main() {
    
    testStdEllipticK();
    testOwnEllipticK();
    
}
