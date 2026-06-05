/****
 * @date Created on 2026-06-05 at 10:16:31 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing implementation of the base class "Cavity".
 ***/
#include "Cavity.hpp"
#include <iostream>

/**
 * Constructor of the Cavity base object.
 */
Cavity::Cavity(const CavityParameters& param) :
        nseed(param.nseed),
        nfreq(param.nfreq),
        freqmin(param.freqmin),
        freqmax(param.freqmax)
    {
    
    // Check for possible errors:
    if (nseed <= 0) {
        std::string msg = "In Cavity(): Invalid number of seeds " + std::to_string(nseed) + ", must be positive.";
        throw std::invalid_argument(msg);
    }
    else if (nfreq <= 0) {
        std::string msg = "In Cavity(): Invalid number of frequencies " + std::to_string(nfreq) + ", must be positive.";
        throw std::invalid_argument(msg);
    }
    else if (freqmin >= freqmax) {
        std::cout << TAG_WARN << "Invalid frequency interval [" << freqmin << ", " << freqmax << "], setting Nfreq=1, and taking Freq=" << (freqmin + freqmax)/2. << " (central value)...";
        nfreq = 1;
    }
    
}

/***********************************
 * GETTERS/SETTERS
 **********************************/

/**
 * Returns the number of random realizations of the disorder over which to average.
 */
int Cavity::getNSeed() const {
    return nseed;
}

/**
 * Returns the number of frequencies over which to average.
 */
int Cavity::getNFreq() const {
    return nfreq;
}

/**
 * Returns the minimum frequency.
 */
double Cavity::getFreqMin() const {
    return freqmin;
}

/**
 * Returns the maximum frequency.
 */
double Cavity::getFreqMax() const {
    return freqmax;
}

/**
 * Returns the frequency at the center of the frequency interval.
 */
double Cavity::getMeanFreq() const {
    return (freqmin + freqmax)/2.;
}

/**
 * Returns the frequency associated with the given index between 0 and nfreq-1.
 */
double Cavity::getFreq(const int ifreq) const {
    if (ifreq < 0 || ifreq >= nfreq) {
        std::string msg = "In getFreq(): Invalid frequency index " + std::to_string(ifreq) + ", expected in 0.." + std::to_string(nfreq-1) + ".";
        throw std::invalid_argument(msg);
    }
    double freq;
    if (nfreq == 1) {
        freq = (freqmin + freqmax)/2.;
    }
    else {
        freq = freqmin + ((freqmax - freqmin)*ifreq)/(nfreq - 1);
    }
    return freq;
}

/**
 * Returns the seed of the realization of the disorder associated with the given index between 0 and nseed-1.
 */
uint64_t Cavity::getSeed(const int iseed) const {
    if (iseed < 0 || iseed >= nseed) {
        std::string msg = "In getSeed(): Invalid seed index " + std::to_string(iseed) + ", expected in 0.." + std::to_string(nseed-1) + ".";
        throw std::invalid_argument(msg);
    }
    return seed0 + iseed;
}

/***********************************
 * PRINTING METHODS
 **********************************/

/**
 * Print the summary to stdout. This method is the same for all cavities.
 */
void Cavity::printSummary() const {
    for (const std::string& line : summary()) {// Loop over the lines of the summary.
        std::cout << TAG_INFO << line << "\n";
    }
}
