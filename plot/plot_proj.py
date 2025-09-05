#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-09-05 at 14:15:10 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot the projection of a 2D scalar field defined over a square lattice on the horizontal direction. The input data must have the form [x, y, f(x, z)].
import sys, os, datetime, csv
import numpy as np
import compile_tikz as ct

def plot_proj(args):
    """
    Plots the component args[1]='column_name' of the given CSV file args[2]='field_file' projected onto the horizontal (x) coordinate and averaged.
    The data in the field file are assumed to have the form [x, y, north, south, east, west, f1, f2, ..., fn], where the components 'x' and 'y' 
    are assumed to be integers, and the cardinal directions are the indices of nearest neighbors or boundary conditions.
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 3):
        print(ct.TAG_ERROR + "Invalid number of arguments, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " COLUMN_NAME FIELD_FILE")
        return 1
    
    column_name = args[1]  ## Interpret arg #1 as the name of the column in the field file.
    field_file = args[2]   ## Interpret arg #2 as the name of the field file.
    file_path = os.path.splitext(field_file)[0]  ## The file path is the filename without its extension (used to write new files). 
    
    try:
        fp = open(field_file, 'r')
    except IOError as e:
        print(ct.TAG_ERROR + "Field file '" + field_file + "' not found, aborting now...")
        return 1
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    data_header = ct.get_header(fp, '%')
    holscat = float(ct.get_value_in_string("h/lscat", data_header))  ## Extract the ratio h/lscat.
    
    ## Extract the bounds of the system from the array of points:
    point = np.asarray([(int(p['x']), int(p['y'])) for p in data], dtype=int)
    xmin = point[:, 0].min()
    xmax = point[:, 0].max()
    ymin = point[:, 1].min()
    ymax = point[:, 1].max()
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1
    matrix = np.full((ny, nx), np.nan)
    
    ## Loop on the points to write in 'matrix':
    for p in data:
        i = ymax - int(p['y'])
        j = int(p['x']) - xmin
        matrix[i, j] = float(p[column_name])
    
    ## Compute the average projection along the 'x' direction:
    proj = np.zeros((nx, 2))
    for ix in range(nx):
        proj[ix, 0] = xmin + (xmax-xmin)*(float(ix)/(nx-1))
    
    proj[:, 1] = np.nanmean(matrix, axis=0)  ## Compute the mean ignoring NaN values.
    
    ## Write the projection data in a string:
    linelen = 3  ## Number of points on each line (arbitrary but not too large).
    proj_string = ""
    for i in range(nx):
        proj_string += "(" + str(proj[i, 0]) + ", " + str(proj[i, 1]) + ") "
        if (i%linelen == linelen-1):
            proj_string += "\n\t"
    
    ## Write the TikZ code and compile the result:
    tikz_code = """%% Generated on {timestamp} by {my_program} {my_copyright}
{data_header}
\\begin{{tikzpicture}}%
\\begin{{axis}}[%
    title={{{title}}},
    xlabel={{{xlabel}}},
    ylabel={{{ylabel}}},
    xmin={xmin}, xmax={xmax},
    xticklabel={{\\pgfmathparse{{{holscat}*\\tick}}$\\pgfmathprintnumber[fixed relative, precision=3]{{\\pgfmathresult}}$}}, %% Rescale ticks to get x/lscat = (x/h) * (h/lscat), with h/lscat={holscat}.
    unbounded coords=jump,  %% Discard NaN's and negative entries.
    clip marker paths=true, %% Clips the marks out of the axis frame.
    clip mode=individual,   %% Ensure the marks do not overlay the other curves.
]%
\\addplot[black, thick, line join=bevel] coordinates {{%% 
\t{proj_string}
}};
\\end{{axis}}%
\\end{{tikzpicture}}%""".format(
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = ct.MY_COPYRIGHT,
        data_header = data_header,
        title = "\\textbf{Cmd:} \\detokenize{"+ " ".join(args) + "}",
        xlabel = "$x/\\ell$",
        ylabel = "\\detokenize{" + column_name + "}",
        holscat = holscat,
        xmin   = xmin,
        xmax   = xmax,
        proj_string = proj_string
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + '_proj.tikz'
    print(ct.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    ct.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    return 0

if (__name__ == '__main__'):
    exit(plot_proj(sys.argv))
