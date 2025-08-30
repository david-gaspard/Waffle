/****
 * @date Created on 2025-08-27 at 13:41:09 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing some tests for the "map" object and other related data structures.
 * Some code snippets are from: https://en.cppreference.com/w/cpp/container/map.html
 ***/
#include "Constants.hpp"
#include "BaseTools.hpp"
#include <vector>
#include <map>
#include <list>
//#include <forward_list>
#include <algorithm>
#include <iostream>

/**************************************************************************************************
 * SPEED TEST FOR SORTED INSERTION IN A MAP
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
 * SPEED TEST FOR SORTED INSERTION IN A VECTOR OF PAIRS
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

/**************************************************************************************************
 * SPEED TEST FOR PUSH_BACK IN A LINKED LIST OF PAIRS + SORT
 **************************************************************************************************/

/**
 * Print a vector of pairs to std output.
 */
void printList(const std::string& name, const std::list<std::pair<uint32_t, uint32_t>>& ls) {
    std::cout << TAG_INFO << name << " = {";
    for (const auto& p : ls) {
        std::cout << " [" << p.first << "]=" << p.second << " ";
    }
    std::cout << "}\n";
}

/**
 * Speed test of the insertion in a linked list of pairs + sort.
 */
double speedInsertListSort(const uint32_t nmax, const bool verbose) {
    std::cout << "====== SPEED TEST : INSERT LIST + SORT ======\n";
    
    std::list<std::pair<uint32_t, uint32_t>> ls;
    uint32_t k = 1;
    
    std::pair<uint32_t, uint32_t> p = {0, 0};
    
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (uint32_t i = 1; i <= nmax; i++) {// Long loop to test insertion time.
        k = (22695477*k + 1) % 2147483648;  // Random key using LCG (see: https://en.wikipedia.org/wiki/Linear_congruential_generator).
        p.first = k;  // Setup the key-value pair.
        p.second = i;
        ls.push_back(p);
    }
    
    double ctime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()/1e6;
    std::cout << TAG_INFO << "Insertion = " << ctime << " s.\n";
    
    ls.sort(compareKey);
    
    ctime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()/1e6 - ctime;
    
    std::cout << TAG_INFO << "Sort = " << ctime << " s.\n";
    
    if (verbose) printList("List", ls);
    
    return endProgressBar(start);
}

/**************************************************************************************************
 * SPEED TEST FOR PUSH_BACK IN A VECTOR OF PAIRS + SORT
 **************************************************************************************************/

/**
 * Speed test of the insertion in a linked list of pairs + sort.
 */
double speedInsertVectorSort(const uint32_t nmax, const bool verbose) {
    std::cout << "====== SPEED TEST : INSERT VECTOR + SORT ======\n";
    
    std::vector<std::pair<uint32_t, uint32_t>> v;
    uint32_t k = 1;
    
    std::pair<uint32_t, uint32_t> p = {0, 0};
    
    const auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    for (uint32_t i = 1; i <= nmax; i++) {// Long loop to test insertion time.
        k = (22695477*k + 1) % 2147483648;  // Random key using LCG (see: https://en.wikipedia.org/wiki/Linear_congruential_generator).
        p.first = k;  // Setup the key-value pair.
        p.second = i;
        v.push_back(p);
    }
    
    double ctime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()/1e6;
    std::cout << TAG_INFO << "Insertion = " << ctime << " s.\n";
    
    //v.sort(compareKey);
    std::sort(v.begin(), v.end(), compareKey);
    
    ctime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()/1e6 - ctime;
    
    std::cout << TAG_INFO << "Sort = " << ctime << " s.\n";
    
    if (verbose) printVector("Vector", v);
    
    return endProgressBar(start);
}

/**
 * Main function of this test.
 */
int main(int argc, char** argv) {
    
    const uint32_t nmax = 1e6;  // Quick lookup: nmax=100. Quick test: nmax=1e5. Big test: nmax=1e6.
    const bool verbose = false;
    
    std::cout << TAG_INFO << "Running speed tests for nmax=" << nmax << "...\n";
    
    const double best_time = speedInsertVectorSort(nmax, verbose);
    
    double ctime;
    ctime = speedInsertListSort(nmax, verbose);
    std::cout << TAG_INFO << ctime/best_time << " times slower than speedInsertVectorSort()...\n";
    
    ctime = speedInsertMap(nmax, verbose);
    std::cout << TAG_INFO << ctime/best_time << " times slower than speedInsertVectorSort()...\n";
    
    ctime = speedInsertVector(nmax, verbose);
    std::cout << TAG_INFO << ctime/best_time << " times slower than speedInsertVectorSort()...\n";
    
    /**
     * TYPICAL OUTPUT FOR NMAX = 10^6 ON A DELL LAPTOP WITH INTEL CORE i7-12800H (EXECUTED IN SERIAL) USING G++ -O2:
     * 
     * [INFO] Running speed tests for nmax=1000000...
     * ====== SPEED TEST : INSERT VECTOR + SORT ======
     * [INFO] Insertion = 0.009964 s.
     * [INFO] Sort = 0.093797 s.
     * [EXEC] Done in 0.103781 s.
     * ====== SPEED TEST : INSERT LIST + SORT ======
     * [INFO] Insertion = 0.026413 s.
     * [INFO] Sort = 0.218295 s.
     * [EXEC] Done in 0.244729 s.
     * [INFO] 2.35813 times slower than speedInsertVectorSort()...
     * ====== SPEED TEST : INSERT MAP ======
     * [EXEC] Done in 0.460005 s.
     * [INFO] 4.43246 times slower than speedInsertVectorSort()...
     * ====== SPEED TEST : INSERT VECTOR ======
     * [EXEC] Done in 199.951 s.
     * [INFO] 1926.66 times slower than speedInsertVectorSort()...
     */
    
    return 0;
}
