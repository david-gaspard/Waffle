/****
 * @date Created on 2025-08-27 at 13:41:09 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing some tests for the "map" object and other related data structures.
 * Some code snippets are from: https://en.cppreference.com/w/cpp/container/map.html
 ***/
#include "Constants.hpp"
#include "BaseTools.hpp"
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

/**************************************************************************************************
 * SPEED TEST FOR INSERTION IN A MAP
 **************************************************************************************************/

/**
 * Print the map to std output.
 */
void printMap(const std::string& name, const std::map<uint32_t, uint32_t>& m) {
    std::cout << TAG_INFO << name << " = {";
    for (const auto& [key, value] : m) {
        std::cout << " [" << key << "]=" << value << " ";
    }
    std::cout << "}\n";
}

/**
 * Speed test of the insertion in "map".
 */
double speedInsertMap(const uint32_t nmax, const bool verbose) {
    std::cout << "====== SPEED TEST : INSERT MAP ======\n";
    
    std::map<uint32_t, uint32_t> m;
    uint32_t k = 1;
    
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (uint32_t i = 1; i <= nmax; i++) {// Long loop to test insertion time.
        k = (22695477*k + 1) % 2147483648;  // Random key using LCG (see: https://en.wikipedia.org/wiki/Linear_congruential_generator).
        m[k] = i;
    }
    
    if (verbose) printMap("Map", m);
    
    return endProgressBar(start);
}

/**************************************************************************************************
 * SPEED TEST FOR INSERTION IN A VECTOR OF PAIRS
 **************************************************************************************************/

/**
 * Print a vector of pairs to std output.
 */
void printVector(const std::string& name, const std::vector<std::pair<uint32_t, uint32_t>>& v) {
    std::cout << TAG_INFO << name << " = {";
    for (const auto& p : v) {
        std::cout << " [" << p.first << "]=" << p.second << " ";
    }
    std::cout << "}\n";
}

/**
 * Comparison function for the vector of pairs used to sort the vector of pairs by keys.
 */
bool compareKey(const std::pair<uint32_t, uint32_t>& pair1, const std::pair<uint32_t, uint32_t>& pair2) {
    return pair1.first < pair2.first;
}

/**
 * Speed test of the insertion in a vector of pairs.
 */
double speedInsertVector(const uint32_t nmax, const bool verbose) {
    std::cout << "====== SPEED TEST : INSERT VECTOR ======\n";
    
    std::vector<std::pair<uint32_t, uint32_t>> v;
    uint32_t k = 1;
    
    std::pair<uint32_t, uint32_t> p = {0, 0};
    
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (uint32_t i = 1; i <= nmax; i++) {// Long loop to test insertion time.
        k = (22695477*k + 1) % 2147483648;  // Random key using LCG (see: https://en.wikipedia.org/wiki/Linear_congruential_generator).
        p.first = k;  // Setup the key-value pair.
        p.second = i;
        v.insert(std::lower_bound(v.begin(), v.end(), p, compareKey), p);
    }
    
    if (verbose) printVector("Vector", v);
    
    return endProgressBar(start);
}

/**
 * Main function of this test.
 */
int main(int argc, char** argv) {
    
    const uint32_t nmax = 100000;
    const bool verbose = false;
    
    speedInsertMap(nmax, verbose);
    speedInsertVector(nmax, verbose);
    
    /**
     * TYPICAL OUTPUT FOR NMAX = 10^6 :
     * 
     * ====== SPEED TEST : INSERT MAP ======
     * [EXEC] Done in 0.401139 s.
     * ====== SPEED TEST : INSERT VECTOR ======
     * [EXEC] Done in 201.424 s.
     * 
     * This means a speedup of ~500 using maps for large lists.
     */
    
    return 0;
}
