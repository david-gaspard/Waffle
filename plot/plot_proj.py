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
    if (len(args) != 4):
        print(ct.TAG_ERROR + "Invalid number of arguments, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " COLUMN_NAME UNIT_LENGTH FIELD_FILE")
        return 1
    
    column_name = args[1]  ## Interpret arg #1 as the name of the column in the field file.
    unit_length = float(eval(args[2]))  ## Interpret arg #6 as the value of h/lscat, or more generally the unit length.
    field_file = args[3]   ## Interpret arg #2 as the name of the field file.
    file_path = os.path.splitext(field_file)[0] + "_" + column_name + "_proj"  ## The file path will be used to write new files.
    
    try:
        fp = open(field_file, 'r')
    except IOError as e:
        print(ct.TAG_ERROR + "Field file '" + field_file + "' not found, aborting now...")
        return 1
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    data_header = ct.get_header(fp, '%')
    
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
        proj[ix, 0] = ((xmax-xmin)*(float(ix)/(nx-1)) + 0.5) * unit_length
    
    proj[:, 1] = np.nanmean(matrix, axis=0)  ## Compute the mean ignoring NaN values.
    
    ## Write the projected data in a CSV file:
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
    np.savetxt(csv_file, proj, fmt='%.16g', delimiter=", ", header=csv_header, comments="")
    
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
        title = "\\textbf{Cmd:} \\detokenize{"+ " ".join(args) + "}",
        xlabel = "$x/\\ell$",
        ylabel = "\\detokenize{" + column_name + "}",
        xmin   = proj[:, 0].min(),
        xmax   = proj[:, 0].max(),
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
    exit(plot_proj(sys.argv))
