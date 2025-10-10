#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-09-16 at 13:18:17 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot a cut through a 2D scalar field defined over a square lattice. The input data must have the form [x, y, f(x, z)].
import sys, os, datetime, csv
import numpy as np
import compile_tikz as ct

def plot_cut_region(args):
    """
    Plots the component args[1]='column_name' of the scalar field defined in the given CSV file args[4]='field_file' averaged horizontally or vertically
    in the rectangular region xmin, xmax, ymin, ymax. Positions are expressed in the coordinates of the mesh points.
    The data in the field file are assumed to have the form [x, y, north, south, east, west, f1, f2, ..., fn], where the components 'x' and 'y' 
    are assumed to be integers, and the cardinal directions are the indices of nearest neighbors or boundary conditions.
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 9):
        print(ct.TAG_ERROR + "Invalid number of arguments, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " COLUMN_NAME XMIN XMAX YMIN YMAX DIRECTION UNIT_LENGTH FIELD_FILE")
        return 1
    
    column_name = args[1]  ## Name of the column in the field file.
    xmin = float(args[2])
    xmax = float(args[3])
    ymin = float(args[4])
    ymax = float(args[5])
    direction = args[6].lower()
    unit_length = float(eval(args[7]))  ## Value of the unit length.
    field_file = args[8]   ## Interpret arg #7 as the name of the field file.
    file_path = os.path.splitext(field_file)[0] + "_" + column_name + "_cut" + direction  ## The file path will be used to write new files.
    
    try:
        fp = open(field_file, 'r')
    except IOError as e:
        print(ct.TAG_ERROR + "Field file '" + field_file + "' not found, aborting now...")
        return 1
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    data_header = ct.get_header(fp, '%')
    
    ## Count the expected number of points in the cut:
    if (direction == "x"):
        npcut = int(np.floor(xmax) - np.ceil(xmin) + 1)
    elif (direction == "y"):
        npcut = int(np.floor(ymax) - np.ceil(ymin) + 1)
    else:
        print(ct.TAG_ERROR + "Invalid direction argument, received '" + direction + "', expected 'x' or 'y', aborting...")
        return 1
    
    cut = np.zeros((npcut, 2))
    nsample = np.zeros(npcut)
    
    ## Loop over the point to compute the cut:
    for p in data:
        x = float(p['x'])
        y = float(p['y'])
        if (xmin <= x and x <= xmax and ymin <= y and y <= ymax):
            if (direction == 'x'):
                i = int(np.floor(x) - np.ceil(xmin))
                cut[i, 0] = x
            else:
                i = int(np.floor(y) - np.ceil(ymin))
                cut[i, 0] = y
            cut[i, 1] += float(p[column_name])
            nsample[i] += 1
    
    ## Normalize the cut abscissa and ordinate:
    cut[:, 0] = (cut[:, 0] - cut[0, 0] + 0.5)*unit_length
    cut[:, 1] /= nsample
    
    ## Write the cut data in a CSV file:
    csv_file = file_path + ".csv"
    csv_header = """%% Generated on {timestamp} by {my_program} {my_copyright}
{data_header}
%% Command: {command}
x, {column_name}""".format(
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = ct.MY_COPYRIGHT,
        data_header = data_header,
        command = " ".join(args),
        column_name = column_name
    )
    np.savetxt(csv_file, cut, fmt='%.16g', delimiter=", ", header=csv_header, comments="")
    
    ## Write the TikZ code and compile the result:
    tikz_code = """%% Generated on {timestamp} by {my_program} {my_copyright}
{data_header}
\\begin{{tikzpicture}}%
\\begin{{axis}}[%
    title={{{title}}},
    xlabel={{{xlabel}}},
    ylabel={{{ylabel}}},
    xmin={xmin}, xmax={xmax},
    ymin=0,
    %%ymode=log,
    unbounded coords=jump,  %% Discard NaN's and negative entries.
    clip marker paths=true, %% Clips the marks out of the axis frame.
    clip mode=individual,   %% Ensure the marks do not overlay the other curves.
]%
\\addplot[black, thick, line join=bevel] table[x=x, y={column_name}]{{{csv_file}}};
\\end{{axis}}%
\\end{{tikzpicture}}%""".format(
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = ct.MY_COPYRIGHT,
        data_header = data_header,
        title  = "\\textbf{Cmd:} \\detokenize{"+ " ".join(args) + "}",
        xlabel = "$x/\\ell$",
        ylabel = "\\detokenize{" + column_name + "}",
        xmin   = cut[0, 0] - 0.5*unit_length,
        xmax   = cut[-1, 0] + 0.5*unit_length,
        column_name = column_name,
        csv_file = "\\jobname.csv"
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + ".tikz"
    print(ct.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    ct.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    return 0


if (__name__ == '__main__'):
    exit(plot_cut_region(sys.argv))
