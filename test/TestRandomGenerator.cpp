/****
 * @date Created on 2025-09-24 at 18:06:54 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the random generator.
 ***/
#include "Constants.hpp"
#include <random>
#include <iostream>

/**
 * Returns the average of the given vector.
 */
double average(const std::vector<double>& sample) {
    
    if (sample.empty()) {
        return 0;
    }
    
    return std::reduce(sample.begin(), sample.end())/static_cast<double>(sample.size());
}

/**
 * Returns the variance of the given vector.
 */
double variance(const std::vector<double>& sample) {
    
    if (sample.empty()) {
        return 0;
    }
    
    const double avg = average(sample);
    double var = 0.;
    
    for (const double s : sample) {
        var += (s-avg)*(s-avg);
    }
    
    return var/static_cast<double>(sample.size());
}

/**
 * Count the number of samples in a given interval.
 */
double count(const std::vector<double>& sample, const double smin, const double smax) {
    
    int cnt = 0;
    
    for (const double s : sample) {
        if (smin < s && s < smax) cnt++;
    }
    
    return cnt;
}

/**
 * Check that the distribution of samples of the normal distribution follows the correct distribution.
 */
void testNormalHisto() {
    std::cout << "====== TEST NORMAL HISTO ======\n";
    const int nsample = 1e6; // Number of samples.
    const uint64_t seed = 1; // Seed used for reproducibility.
    const double mu = 0.;     // Target mean.
    const double sigma = 3.;  // Target standard deviation.
    
    std::mt19937_64 rng;  // Instantiate the standard Mersenne Twister random number generator (64-bit-return version).
    rng.seed(seed);       // Initialize the random generator with the given seed.
    std::normal_distribution<double> random_normal(mu, sigma);
    
    std::vector<double> sample(nsample, 0.);
    
    for (int i = 0; i < nsample; i++) {//Loop over the points of the mesh.
        sample.at(i) = random_normal(rng);
    }
    
    const double avg = average(sample);
    
    const double var = variance(sample);
    const double var_expc = sigma*sigma;
    
    const int cnt_p = count(sample, mu, mu+sigma);
    const int cnt_m = count(sample, mu-sigma, mu);
    const double cnt_expc = 0.341345*nsample;
    
    std::cout << TAG_INFO << "Average  = " << avg << " (expected " << mu << ").\n";
    std::cout << TAG_INFO << "Variance = " << var << " (expected " << var_expc << ").\n";
    std::cout << TAG_INFO << "Count in [mu, mu+sigma] = " << cnt_p << " (expected " << cnt_expc << "), error = " << 100.*std::abs(cnt_p - cnt_expc)/cnt_expc << "%.\n";
    std::cout << TAG_INFO << "Count in [mu-sigma, mu] = " << cnt_m << " (expected " << cnt_expc << "), error = " << 100.*std::abs(cnt_m - cnt_expc)/cnt_expc << "%.\n";
}

/**
 * Main function of this test.
 */
int main(int argc, char** argv) {
    
    testNormalHisto();
    
    return 0;
}
