/****
 * @date Created on 2025-07-09 at 16:09:04 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing basic tools for file manipulation.
 ***/
#include "BaseTools.hpp"
#include "Constants.hpp"
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <iostream>

/**
 * Write a one-liner timestamp in the given stream with the given prefix (typically a comment character when writing to a file).
 */
void writeTimestamp(std::ofstream& ofs, const char* prefix) {
    std::time_t time_now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ofs << prefix << "Computed on " << std::put_time(std::localtime(&time_now), "%F at %T %z") << " by " << PROGRAM_COPYRIGHT << "\n";
}

/**
 * Create a unique filename at the given "path" and with the given "suffix". The unique filename found is written in the string "unique_filename".
 * The format is: <path><i><suffix> where <i> is an integer at least equal to 1 which is incremented until the filename is unique.
 * This function ensures data is not overwritten. If no unique filename is found, then prints a warning message.
 * This functions also creates the appropriate subdirectories if necessary, and warns the user if so.
 */
void uniqueFilename(const std::string& path, const std::string& suffix, std::string& unique_filename) {
    
    // 1. Ensure the output directory exists:
    std::string outputdir(path.substr(0, path.find_last_of("/")));
    std::filesystem::create_directories(outputdir); // Create the directory if necessary.
    
    // 2. Create a unique filename:
    for (int i = 1; i <= MAX_FILENO; i++) {
        unique_filename = path + std::to_string(i) + suffix;
        if (!std::filesystem::exists(unique_filename)) return; // If the file does not exist, then exits the function.
    }
    
    std::cout << TAG_WARN << "Could not find unique filename: '" << unique_filename << "'. Data will be lost...\n";
}

/**
 * Returns a string version of the given time duration in seconds in the format HH:MM:SS.
 */
const std::string timeToString(double time) {
    
    int time_hours   = (int) std::floor(time/3600.);
    time -= 3600.*time_hours;
    int time_minutes = (int) std::floor(time/60.);
    time -= 60.*time_minutes;
    int time_seconds = (int) std::floor(time);
    
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << time_hours << ":" << std::setw(2) << time_minutes << ":" << std::setw(2) << time_seconds;
    
    return ss.str();
}

/**
 * Print a progress bar indicator to std output.
 * 
 * Arguments:
 * 
 * cjob  = Number of completed jobs.
 * njob  = Total number of jobs to be completed.
 * msg   = Message to be shown in std output. Typically the name of the task.
 * start = Time reference to some fixed point in the past.
 */
void printProgressBar(const int64_t cjob, const int64_t njob, const std::string& msg, const std::chrono::steady_clock::time_point& start) {
    
    double elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()/1e6;
    
    // 1. Compute the estimated remaining time, expected time of arrival (ETA):
    double eta = (elapsed_time * (njob - cjob))/cjob;  // ETA in number of seconds.
    std::string eta_str = timeToString(eta);
    
    // 2. Compute the progress bar:
    int len = (int) std::ceil(((double) BAR_LENGTH * cjob)/njob); // Length of the filled part of the progress bar.
    int percent = (int) std::ceil((100.*cjob)/njob); // Completed jobs in percent.
    int njobw = (int) std::ceil(std::log10(njob)); // Width of the number of jobs (in base-10 decimals).
    
    std::cout << std::setfill(' ')
        << TAG_EXEC << msg << " [" << std::string(len, '#') << std::string(BAR_LENGTH - len, ' ') << "] " 
        << std::setw(3) << percent << "% (" << std::setw(njobw) << cjob << "/" << std::setw(njobw) << njob << ")  ETA " << eta_str << "  \r" << std::flush;
    
}

/**
 * Finalize the progress bar by computing and returning the total elapsed time (in seconds).
 */
double endProgressBar(const std::chrono::steady_clock::time_point& start) {
    
    double elapsed_sec = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()/1e6;
    
    std::cout << "\n" << TAG_EXEC << "Done in " << elapsed_sec << " s.\n";
    
    return elapsed_sec;
}

/**
 * Convert the given "value" to string using given precision "prec".
 */
std::string to_string_prec(const double value, const int prec) {
    std::stringstream ss;
    ss.precision(prec);
    ss << value;
    return ss.str();
}
