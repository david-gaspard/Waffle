/****
 * @date Created on 2025-08-12 at 11:07:41 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing raster graphics utilities, especially reading/writing raster images.
 ***/
#include "Image.hpp"
#include <png.h>  // Official libpng library: http://www.libpng.org/pub/png/libpng.html
#include <filesystem>
#include <iostream>

/**************************************************************************************************
 * RASTER IMAGE MANIPULATION
 **************************************************************************************************/

/**
 * Constructor of the Image object. Note that it is (height, width) and NOT (width, height) to be consistent with (i, j) indices in matrix conventions.
 */
Image::Image(const int height, const int width) {
    if (height <= 0) {
        throw std::invalid_argument("In Image(): Height cannot be negative.");
    }
    else if (width <= 0) {
        throw std::invalid_argument("In Image(): Width cannot be negative.");
    }
    //std::cout << TAG_INFO << "Creating Image...\n";
    this->width = width;
    this->height = height;
    data = new Color[width*height](); // Initialize to zero (black image).
}

/**
 * Constructor of the Image object from a file.
 */
Image::Image(const std::string& filename) {
    //std::cout << TAG_INFO << "Creating Image from file '" << filename << "'..." << std::endl;
    
    if (not std::filesystem::exists(filename)) throw std::invalid_argument("In Image(): File '" + filename + "' not found");
    
    FILE *fp = fopen(filename.c_str(), "rb");
    
    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) throw std::runtime_error("In Image(): Aborted because !png");
    
    png_infop info = png_create_info_struct(png);
    if (!info) throw std::runtime_error("In Image(): Aborted because !info");
    if (setjmp(png_jmpbuf(png))) throw std::runtime_error("In Image(): Aborted because setjmp()");
    
    png_init_io(png, fp);
    png_read_info(png, info);
    width  = png_get_image_width(png, info);
    height = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth  = png_get_bit_depth(png, info);
    
    if (bit_depth == 16)  // Convert 16bit to 8bit depth.
        png_set_strip_16(png);
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png);
    if(color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png);
    if (color_type == PNG_COLOR_TYPE_PALETTE) 
        png_set_palette_to_rgb(png);
    if (color_type & PNG_COLOR_MASK_ALPHA)  // Strip alpha channel.
        png_set_strip_alpha(png);
    
    png_read_update_info(png, info);
    
    const int nchannels = png_get_channels(png, info); // Get the number of channels...
    if (nchannels != 3) {// Check that the number of channels is 3 as expected.
        std::string msg = "In Image(): Number of channels must be 3 (received " + std::to_string(nchannels) + ")...";
        throw std::invalid_argument(msg);
    }
    auto row_pointers = new png_bytep[height];
    for (int i = 0; i < height; i++) {// Loop over the rows.
        row_pointers[i] = new png_byte[nchannels*width];
    }
    
    png_read_image(png, row_pointers); // Read the PNG file. Data is now in "row_pointers".
    fclose(fp);
    
    // Now convert the data in "row_pointers" to the RGB data of the present image:
    data = new Color[width*height]();
    for (int i = 0; i < height; i++) {// Loop over the rows.
        for (int j = 0; j < width; j++) {// Loop over the pixels of the row.
            data[i*width + j] = {
                row_pointers[i][nchannels*j], 
                row_pointers[i][nchannels*j+1], 
                row_pointers[i][nchannels*j+2]
            };
        }
        delete[] row_pointers[i];
    }
    delete[] row_pointers;
}

/**
 * Destructor of the Image object.
 */
Image::~Image() {
    //std::cout << TAG_INFO << "Deleting Image...\n";
    delete[] data;
}

/**************************************************************************************************
 * GETTERS/SETTERS
 **************************************************************************************************/

/**
 * Returns the height of the image.
 */
int Image::getHeight() const {
    return height;
}

/**
 * Returns the width of the image.
 */
int Image::getWidth() const {
    return width;
}

/**
 * returns the color at position (i, j) in the image.
 * Index "i" is the row index, and "j" the column index, just as in a matrix.
 */
Color& Image::operator()(const int i, const int j) const {
    if (i >= height || j >= width || i < 0 || j < 0) {
        std::string msg = "In Image(): Invalid index (i=" + std::to_string(i) + ", j=" + std::to_string(j) 
                        + ") given height=" + std::to_string(height) + " and width=" + std::to_string(width) + ".";
        throw std::out_of_range(msg);
    }
    return data[i*width + j];
}

/**************************************************************************************************
 * IMAGE EXPORT
 **************************************************************************************************/

/**
 * Save the present image in Netpbm portable pixmap format (PPM).
 */
void Image::savePPM(const std::string& filename) const {
    
    // Open a file in binary mode:
    std::ofstream ofs(filename, std::ios::binary);
    ofs << "P6 " << (uint32_t)width << ' ' << (uint32_t)height << ' ' << (uint32_t)MAX_COLOR << ' ';
    
    for (int l = 0; l < width*height; l++) {// Loop over all pixels.
        ofs << data[l];
    }
    ofs.close(); // Make sure the stream is closed (it should be).
}

/**
 * Save the present image in Portable Network Graphics (PNG) format.
 */
void Image::savePNG(const std::string& filename) const {
    
    FILE* fp = fopen(filename.c_str(), "wb"); // Open the file.
    if (!fp) std::abort();
    
    // Initialize the PNG instance:
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) std::abort();
    
    png_infop info = png_create_info_struct(png);
    if (!info) std::abort();
    if (setjmp(png_jmpbuf(png))) std::abort();
    
    png_init_io(png, fp);
    
    png_set_IHDR(
        png,
        info,
        width, height,
        8,
        PNG_COLOR_TYPE_RGB,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);
    
    // Convert the data to PNG row_pointers:
    const int nchannels = 3; // There are always 3 channels (red, green, blue).
    Color color;
    auto row_pointers = new png_bytep[height];
    for (int i = 0; i < height; i++) {// Loop over the rows.
        row_pointers[i] = new png_byte[nchannels*width];
        for (int j = 0; j < width; j++) {// Loop over the pixels of the row.
            color = data[i*width + j];
            row_pointers[i][nchannels*j]   = color.r;
            row_pointers[i][nchannels*j+1] = color.g;
            row_pointers[i][nchannels*j+2] = color.b;
        }
    }
    
    // Write the PNG file:
    png_write_image(png, row_pointers);
    
    // End the write process:
    png_write_end(png, NULL);
    for (int i = 0; i < height; i++) {
        delete[] row_pointers[i];
    }
    delete[] row_pointers;
    fclose(fp);
}

/**************************************************************************************************
 * IMAGE TRANSFORMATIONS
 **************************************************************************************************/

/**
 * Fill the present image with the given color.
 */
void Image::fill(const Color color) {
    std::fill(data, data + width*height, color);
}

/**
 * Invert the colors of the present image.
 */
void Image::invert() {
    Color color;
    for (int l = 0; l < width*height; l++) {// Loop over the pixels.
        color = data[l];
        data[l] = {
            static_cast<uint8_t>(MAX_COLOR - color.r), 
            static_cast<uint8_t>(MAX_COLOR - color.g), 
            static_cast<uint8_t>(MAX_COLOR - color.b)
        };
    }
}
