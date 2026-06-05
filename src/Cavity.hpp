/****
 * @date Created on 2026-06-04 at 17:37:26 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing header for the base "Cavity" class which contains the parameters defining any cavity geometry.
 ***/
#ifndef _CAVITY_H
#define _CAVITY_H
#include "WaveSystem.hpp"

/**
 * Structure containing the cavity parameters to pass to the constructor.
 */
struct CavityParameters {
    
    int nseed;       // Number of realizations of the disorder. Recommended is 1000.
    int nfreq;       // Number of frequencies in the frequency range. The interval [freqmin, freqmax] is covered uniformly at constant steps. Nfreq=1 selects the center of the interval.
    double freqmin;  // Minimum field frequency of the frequency range (in Hz). Default=12.0e9 Hz.
    double freqmax;  // Minimum field frequency of the frequency range (in Hz). Default=15.0e9 Hz.
    
};

/**
 * Base class defining a generic "Cavity" object.
 */
class Cavity {
    
    private:
    
    // Members shared by any cavity:
    const uint64_t seed0 = 1;  // Initial seed of the disorder generation. Always equal to 1.
    int nseed;       // Number of realizations of the disorder. Recommended is 1000.
    int nfreq;       // Number of frequencies in the frequency range. The interval [freqmin, freqmax] is covered uniformly at constant steps. Nfreq=1 selects the center of the interval.
    double freqmin;  // Minimum field frequency of the frequency range (in Hz). Default=12.0e9 Hz.
    double freqmax;  // Minimum field frequency of the frequency range (in Hz). Default=15.0e9 Hz.
    
    public:
    
    // Constructor:
    Cavity(const CavityParameters& param);
    
    // Abstract methods overriden by any cavity:
    virtual std::vector<std::string> summary() const = 0;
    virtual WaveSystem generateSystem(const uint64_t seed, const double freq) const = 0;
    
    // Method which are the same for every cavity:
    int getNSeed() const;
    int getNFreq() const;
    double getFreqMin() const;
    double getFreqMax() const;
    
    double getMeanFreq() const;
    double getFreq(const int ifreq) const;
    uint64_t getSeed(const int iseed) const;
    
    void printSummary() const;
    
};

#endif
