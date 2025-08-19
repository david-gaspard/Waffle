/****
 * @date Created on 2025-08-09 at 20:29:50 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the modulo.
 ***/
#include <iostream>
#include <cmath>

/**
 * Test the modulo for a given number.
 */
void testModulo(const double x) {
    std::cout << "====== TEST MODULO FOR: " << x << " ======\n";
    
    double integer_part;
    double fractional_part = std::modf(x, &integer_part);
    
    std::cout << "[INFO] Integer part:    " << integer_part << "\n";
    std::cout << "[INFO] Fractional part: " << fractional_part << "\n";
}

/**
 * Main function of the test.
 */
int main() {
    
    testModulo(-3.5678);
    testModulo(-2.5678);
    testModulo(-1.5678);
    testModulo( 0.5678);
    testModulo( 1.5678);
    testModulo( 2.5678);
    
    return 0;
}
