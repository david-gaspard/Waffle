/****
 * @date Created on 2025-09-05 at 17:43:53 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the Image class.
 ***/
#include "Image.hpp"
#include <iostream>

/**
 * Create a small images for tests.
 */
Image createImageSmall() {
    Image img(2, 3);
    
    img(0, 0) = COLOR_RED;
    img(0, 1) = COLOR_GREEN;
    img(0, 2) = COLOR_BLUE;
    img(1, 0) = COLOR_YELLOW;
    img(1, 1) = COLOR_WHITE;
    img(1, 2) = COLOR_BLACK;
    
    return img;
}

/**
 * Test the export in Netpbm portable pixmap format (PPM).
 */
void testSavePPM(const Image& img) {
    std::cout << "====== TEST SAVE PPM ======\n";
    const std::string filename = "old/test/image/small.ppm";
    std::cout << TAG_INFO << "Saving image to '" << filename << "'...\n";
    img.savePPM(filename);
}

/**
 * Test the export in Portable Network Graphics (PNG) format.
 */
void testSavePNG(const Image& img) {
    std::cout << "====== TEST SAVE PNG ======\n";
    const std::string filename = "old/test/image/small.png";
    std::cout << TAG_INFO << "Saving image to '" << filename << "'...\n";
    img.savePNG(filename);
}

/**
 * Test opening a PNG file.
 */
void testOpenInvertPNG() {
    std::cout << "====== TEST OPEN/INVERT PNG ======\n";
    
    Image img("old/test/image/small_alpha.png");  // Open the PNG image.
    
    img.invert();  // Invert the colors.
    
    const std::string outfile = "old/test/image/small_invert.png";
    std::cout << TAG_INFO << "Saving image to '" << outfile << "'...\n";
    img.savePNG(outfile);
}

/**
 * Main function of the test of the image class.
 */
int main(int argc, char** argv) {
    
    Image img = createImageSmall();
    testSavePPM(img);
    testSavePNG(img);
    testOpenInvertPNG();
    
    return 0;
}
