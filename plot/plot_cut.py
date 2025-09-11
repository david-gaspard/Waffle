#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-09-04 at 18:53:19 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot a cut through a 2D scalar field defined over a square lattice. The input data must have the form [x, y, f(x, z)].
import sys, os, datetime, csv
import numpy as np
import compile_tikz as ct

def comparePoint(p1, p2):
    """
    Comparison function.
    Return True if point "p1" is located before "p2" assuming the points are stored in column-major ordering.
    """
    p1x = int(p1['x'])
    p1y = int(p1['y'])
    p2x = int(p2['x'])
    p2y = int(p2['y'])
    return p1x < p2x or (p1x == p2x and p1y > p2y)

def find_binary_search(data, x, y):
    """
    Find the point at position (x, y) (couple of integers) in the "data" and return the index.
    This function uses binary search assuming that the points in "data" are sorted in column-major ordering
    (i.e., they satisfy the comparePoint() function above).
    """
    ileft = 0
    iright = len(data) - 1
    ptarget = {'x': str(x), 'y': str(y)}
    
    while ileft <= iright:
        imid = (ileft + iright)//2
        pmid = data[imid]
        
        if not comparePoint(pmid, ptarget) and not comparePoint(ptarget, pmid):
            return imid
        
        if comparePoint(pmid, ptarget):
            ileft = imid + 1
        else:
            iright = imid - 1
    
    return len(data)

def interpolate(data, column_name, x, y):
    """
    Returns the bilinear interpolated value of the given "data" list for the scalar field "column_name" at position (x, y).
    """
    npoint = len(data)  ## Total number of points of the mesh.
    
    ## 1. Determine the position of the theoretical point in the lower left corner:
    x00 = int(np.floor(x))
    y00 = int(np.floor(y))
    
    ## 2. Find the corresponding point in the data using binary search:
    i00 = find_binary_search(data, x00, y00)
    
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
    if (len(args) != 8):
        print(ct.TAG_ERROR + "Invalid number of arguments, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " COLUMN_NAME AX AY BX BY UNIT_LENGTH FIELD_FILE")
        return 1
    
    column_name = args[1]  ## Interpret arg #1 as the name of the column in the field file.
    a = np.asarray((args[2], args[3]), dtype=float) ## Interpret the following arguments as point "a".
    b = np.asarray((args[4], args[5]), dtype=float) ## Interpret the following arguments as point "b".
    unit_length = float(eval(args[6]))  ## Interpret arg #6 as the value of h/lscat, or more generally the unit length.
    field_file = args[7]   ## Interpret arg #7 as the name of the field file.
    file_path = os.path.splitext(field_file)[0] + "_" + column_name + "_cut"  ## The file path will be used to write new files.
    
    try:
        fp = open(field_file, 'r')
    except IOError as e:
        print(ct.TAG_ERROR + "Field file '" + field_file + "' not found, aborting now...")
        return 1
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    data_header = ct.get_header(fp, '%')
    
    ## Construct the cross-sectional cut using bilinear interpolation:
    L = np.linalg.norm(a - b)  ## Length of the path.
    nsub = int(np.round(L))  ## Appropriate number of subdivisions (intervals) to get unit approximately distance between successive points.
    cut = np.zeros((nsub+1, 2))  ## Initialize the data matrix (s*L, f).
    
    for i in range(nsub+1):
        s = float(i)/nsub ## Linear parameter in [0, 1].
        r = a*(1-s) + s*b
        f = interpolate(data, column_name, r[0], r[1])
        cut[i, 0] = s*L*unit_length
        cut[i, 1] = f
    
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
        xmin   = cut[0, 0],
        xmax   = cut[nsub, 0],
        column_name = column_name,
        csv_file = "\\jobname.csv"
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + ".tikz"
    print(ct.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    ct.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    ## TODO: If time permits, maybe also show the corresponding cut from above using plot_map.py...
    return 0


if (__name__ == '__main__'):
    exit(plot_cut(sys.argv))
