#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-09-04 at 18:53:19 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot a cut through a 2D scalar field defined over a square lattice. The input data must have the form [x, y, f(x, z)].
import sys, os, re, datetime, csv
import numpy as np
import matplotlib.pyplot as mplt
import matplotlib.colors as mcol
import compile_tikz
#import plot_map


def interpolate(data, column_name, x, y):
    """
    Returns the bilinear interpolated value of the given "data" list for the scalar field "column_name" at position (x, y).
    """
    npoint = len(data)  ## Total number of points of the mesh.
    
    ## 1. Find the theoretical point in the lower left corner:
    x00 = int(np.floor(x))
    y00 = int(np.floor(y))
    
    ## 2. Find the corresponding values in the data (NB: this operation is O(N), but optimizable to O(log(N)) using binary search):
    i00 = 0
    while (i00 != npoint and not (int(data[i00]['x']) == x00 and int(data[i00]['y']) == y00)):
        i00 += 1
    
    if (i00 == npoint):
        return 0. ## If the point is not found, then returns zero.
    
    f00 = float(data[i00][column_name])
    
    ## 3. Find the three other points:
    i00north = data[i00]['north']
    if (i00north.isdigit()):
        i01 = int(i00north)
        x01 = int(data[i01]['x'])
        y01 = int(data[i01]['y'])
        f01 = float(data[i01][column_name])
    else:
        return 0.
    
    i00east = data[i00]['east']
    if (i00east.isdigit()):
        i10 = int(i00east)
        x10 = int(data[i10]['x'])
        y10 = int(data[i10]['y'])
        f10 = float(data[i10][column_name])
    else:
        return 0.
    
    i10north = data[i10]['north']
    if (i10north.isdigit()):
        i11 = int(i10north)
        x11 = int(data[i11]['x'])
        y11 = int(data[i11]['y'])
        f11 = float(data[i11][column_name])
    else:
        return 0.
    
    ## 4. Use bilinear interpolation:
    xu = x - x00
    yu = y - y00
    
    return f00 * (1-xu) * (1-yu) + f01 * (1-xu) * yu + f10 * xu * (1-yu) + f11 * xu * yu

def plot_cut(args):
    """
    Plots the component args[1]='column_name' along a straight line between positions (ax, ay)  and (bx, by) of the scalar field defined
    in the given CSV file args[4]='field_file'. Positions are expressed in the coordinates of the mesh points.
    The data in the field file are assumed to have the form [x, y, north, south, east, west, f1, f2, ..., fn], where the components 'x' and 'y' 
    are assumed to be integers, and the cardinal directions are the indices of nearest neighbors or boundary conditions.
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 7):
        print(compile_tikz.TAG_ERROR + "Invalid number of arguments, doing nothing...")
        print(compile_tikz.TAG_USAGE + args[0] + " COLUMN_NAME AX AY BX BY FIELD_FILE")
        return 1
    
    column_name = args[1]  ## Interpret arg #1 as the name of the column in the field file.
    a = np.asarray((args[2], args[3]), dtype=float) ## Interpret the following arguments as point "a".
    b = np.asarray((args[4], args[5]), dtype=float) ## Interpret the following arguments as point "b".
    field_file = args[6]   ## Interpret arg #6 as the name of the field file.
    file_path = os.path.splitext(field_file)[0]  ## The file path is the filename without its extension (used to write new files). 
    
    try:
        fp = open(field_file, 'r')
    except IOError as e:
        print(compile_tikz.TAG_ERROR + "Field file '" + field_file + "' not found, aborting now...")
        return 1
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    
    ## Import the header with essential information on the simulation:
    with open(field_file, 'r') as f:
        data_header = "".join([l for l in f if l.startswith("%")]).strip()
    
    ## Construct the cross-sectional cut using bilinear interpolation:
    L = np.linalg.norm(a - b)  ## Length of the path.
    nsub = int(np.round(L))  ## Appropriate number of subdivisions (intervals) to get unit approximately distance between successive points.
    cut = np.zeros((nsub+1, 2))  ## Initialize the data matrix (s*L, f).
    
    for i in range(nsub+1):
        s = float(i)/nsub ## Linear parameter in [0, 1].
        r = a*(1-s) + s*b
        f = interpolate(data, column_name, r[0], r[1])
        cut[i, 0] = s*L
        cut[i, 1] = f
    
    ## Write the cut data in a string:
    cut_string = ""
    for i in range(nsub+1):
        cut_string += "(" + str(cut[i, 0]) + ", " + str(cut[i, 1]) + ") "
    
    ## Write the TikZ code and compile the result:
    tikz_code = """%% Generated on {timestamp} by {my_program} {my_copyright}
{data_header}
\\begin{{tikzpicture}}%
\\begin{{axis}}[%
    title={{\\detokenize{{{data_file}}}}},
    xlabel={{{xlabel}}},
    ylabel={{\\detokenize{{{ylabel}}}}},
    xmin={xmin}, xmax={xmax},
    unbounded coords=jump,  %% Discard NaN's and negative entries.
    clip marker paths=true, %% Clips the marks out of the axis frame.
    clip mode=individual,   %% Ensure the marks do not overlay the other curves.
]%
\\addplot[black, thick] coordinates {{%% 
    {cut_string}
}};
\\end{{axis}}%
\\end{{tikzpicture}}%""".format(
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = compile_tikz.MY_COPYRIGHT,
        data_header = data_header,
        data_file = field_file,
        xlabel = "$x$",
        ylabel = column_name,
        xmin   = cut[0, 0],
        xmax   = cut[nsub, 0],
        cut_string = cut_string
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + '_cut.tikz'
    print(compile_tikz.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    compile_tikz.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    ## TODO: If time permits, maybe also show the corresponding cut from above using plot_map.py...
    
    return 0


if (__name__ == '__main__'):
    exit(plot_cut(sys.argv))
