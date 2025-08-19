/****
 * @date Created on 2025-07-09 at 16:11:42 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing basic tools for file manipulation.
 ***/
#ifndef _BASE_TOOLS_H
#define _BASE_TOOLS_H
#include <fstream>
#include <chrono>

static const int MAX_FILENO = 10000; // Set a maximum file number (safety limit).
static const int BAR_LENGTH = 50;    // Length of the progress bar in number of characters.

void writeTimestamp(std::ofstream& ofs, const char* prefix);
void uniqueFilename(const std::string& path, const std::string& suffix, std::string& unique_filename);
void printProgressBar(const int64_t cjob, const int64_t njob, const std::string& msg, const std::chrono::steady_clock::time_point& start);
double endProgressBar(const std::chrono::steady_clock::time_point& start);
std::string to_string_prec(const double value, const int prec);

#endif
