#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-08-28 at 14:12:25 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot the histogram of the sample data stored in a CSV file.
import sys, os, datetime
import numpy as np
import compile_tikz as ct

def plot_histo(args):
    """
    Plot the histogram of the sample data stored in the given file.
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 2):
        print(ct.TAG_ERROR + "No input file, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " SAMPLE_FILE")
        return 1
    
    sample_file = args[1]
    file_path = os.path.splitext(sample_file)[0]  ## The file path is the filename without its extension (used to write new files). 
    
    ## Import the data:
    #start = time.time()
    samples = np.genfromtxt(sample_file, delimiter=",", comments="%")
    samples = samples[1:, :]
    #end = time.time()
    #print(ct.TAG_INFO + "File imported in", end - start, "s.")
    
    #print(samples[0:6, 0:6])
    nsample = samples.size
    #print(ct.TAG_INFO + "Found", nsample, "samples...")
    
    ## Import the header with essential information on the simulation:
    with open(sample_file, 'r') as fp:
        data_header = ct.get_header(fp, '%')
    
    ## Set up the histogram:
    xmin = 0.  ## Lower bound on the samples.
    xmax = 1.  ## Upper bound on the samples.
    nbin = int(np.ceil(2. * nsample**(1./3)))  ## Use Rice's rule to estimate the number of bins (slightly overestimate the optimal bin number).
    
    ## Use Chebyshev nodes for the bins:
    funbin = lambda i : xmin + (xmax - xmin)*(1. - np.cos(np.pi*i/nbin))/2  ## Bin function for i=[0, 1, ..., nbin]
    binlist = np.array([funbin(i) for i in range(nbin+1)])
    
    ## Compute the histogram:
    histo = np.histogram(samples, bins=binlist, density=True)  ## The sample array is implicitly flattened.
    tavg = np.mean(samples)  ## Compute the average transmission probability (to adjust the bimodal law).
    
    ## Write the histogram data string:
    linelen = 3  ## Number of points on each line (arbitrary but not too large).
    histo_string = ""
    for i in range(histo[0].size):
        histo_string += "(" + str(funbin(i+0.5)) + ", " + str(histo[0][i]) + ") "
        if (i%linelen == linelen-1):
            histo_string += "\n\t"
    
    ## Write the TikZ code:
    tikz_code = """%% Generated on {timestamp} by {my_program} {my_copyright}
{data_header}
\\begin{{tikzpicture}}%
\\pgfmathsetmacro\\tavg{{{tavg}}}%% Average transmission, Tavg ~ 1/[1 + (2/pi) * (L/lscat)]
\\begin{{axis}}[%
    title={{{title}}},
    xlabel={{{xlabel}}},
    ylabel={{{ylabel}}},
    xmin={xmin}, xmax={xmax},
    ymin=0.02, ymax=200,
    ymode=log,
    unbounded coords=jump,  %% Discard NaN's and negative entries.
    clip marker paths=true, %% Clips the marks out of the axis frame.
    clip mode=individual,   %% Ensure the marks do not overlay the other curves.
]%
\\addplot[black!30, smooth, domain=0.001:0.999, samples=64] ({{sin(90*\\x)^2}}, {{\\tavg/(2 * sin(90*\\x)^2 * cos(90*\\x))}}); %% Plot the bimodal distribution.
\\addplot[mark=*, only marks, mark size=0.8] coordinates {{%% 
\t{histo_string}
}};
\\end{{axis}}%
\\end{{tikzpicture}}%""".format(
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = ct.MY_COPYRIGHT,
        data_header = data_header,
        title = "\\textbf{Cmd:} \\detokenize{"+ " ".join(args) + "}",
        xlabel = "Transmission eigenvalue $T$",
        ylabel = "Distribution $\\rho(T)$",
        xmin   = 0,
        xmax   = 1,
        tavg = tavg,
        histo_string = histo_string
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + '.tikz'
    print(ct.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    ct.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    return 0

if (__name__ == '__main__'):
    exit(plot_histo(sys.argv))
