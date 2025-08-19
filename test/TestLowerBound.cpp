/**
 * @date Created on 2025-08-05 at 15:58:49 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * Based on the example from: https://cplusplus.com/reference/algorithm/lower_bound/
 */
#include <iostream>     // std::cout
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector

/**
 * Standard test for lower_bound() and upper_bound().
 * Example from: https://cplusplus.com/reference/algorithm/lower_bound/
 */
int testLowerBound1() {
    std::cout << "====== TEST LOWER BOUND #1 ======\n";
    
    int myints[] = {10,20,30,30,20,10,10,20};
    std::vector<int> v(myints,myints+8);            // 10 20 30 30 20 10 10 20
    std::sort(v.begin(), v.end());                  // 10 10 10 20 20 20 30 30
    
    std::vector<int>::iterator low, up;
    low = std::lower_bound(v.begin(), v.end(), 20); //          ^
    up  = std::upper_bound(v.begin(), v.end(), 20); //                   ^
    
    std::cout << "lower_bound at position " << (low - v.begin()) << ", element " << *low << "\n";
    std::cout << "upper_bound at position " << (up  - v.begin()) << ", element " << *up  << "\n";
    
    return 0;
}

/**
 * Second test for lower_bound() and upper_bound().
 */
int testLowerBound2() {
    std::cout << "====== TEST LOWER BOUND #2 ======\n";
    
    int myints[] = {10,20,30,30,20,10,10,20};
    std::vector<int> v(myints,myints+8);            // 10 20 30 30 20 10 10 20
    
    std::sort(v.begin(), v.end());                  // 10 10 10 20 20 20 30 30
    
    std::vector<int>::iterator low, up;
    low = std::lower_bound(v.begin(), v.end(), 15); //          ^ 
    up  = std::upper_bound(v.begin(), v.end(), 15); //          ^ 
    
    std::cout << "lower_bound at position " << (low - v.begin()) << ", element " << *low << "\n";
    std::cout << "upper_bound at position " << (up  - v.begin()) << ", element " << *up  << "\n";
    
    return 0;
}

/**
 * Third test for lower_bound() and upper_bound().
 */
int testLowerBound3() {
    std::cout << "====== TEST LOWER BOUND #3 ======\n";
    
    std::vector<int> v = {10, 11, 12};
    
    std::sort(v.begin(), v.end());
    
    std::vector<int>::iterator low, up;
    low = std::lower_bound(v.begin(), v.end(), 15);
    up  = std::upper_bound(v.begin(), v.end(), 15);
    
   std::cout << "lower_bound at position " << (low - v.begin()) << ", element " << *low << "\n";
   std::cout << "upper_bound at position " << (up  - v.begin()) << ", element " << *up  << "\n";
    
    auto p = v.insert(low, 15);
    
    std::cout << "v = [";
    for (unsigned int i = 0; i < v.size(); i++) {
        std::cout << " " << v.at(i) << " ";
    }
    std::cout << "]\n";
    
    std::cout << "inserted " << *p << " at position " << p - v.begin() << "\n";
    
    return 0;
}

/**
 * Main function of the test for lower_bound().
 */
int main () {
    
    testLowerBound1();
    testLowerBound2();
    testLowerBound3();
    
    return 0;
}
