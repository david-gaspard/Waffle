#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-07-17 at 12:48:24 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot a scalar field defined over a square lattice with holes. The input data must have the form [x, y, f(x, y)].
import sys, os, datetime, csv
import numpy as np
import matplotlib.pyplot as mplt
import matplotlib.colors as mcol
import compile_tikz as ct
import plot_mesh

## Create the 'sunset' colormap, a luminance-linearized colormap originally created by David Gaspard in June 2024:
SUNSET_COLORS = [[1.00000, 1.00000, 1.00000],
                 [1.00000, 0.99337, 0.65453],
                 [1.00000, 0.96893, 0.46622],
                 [1.00000, 0.93848, 0.36430],
                 [1.00000, 0.90595, 0.29694],
                 [1.00000, 0.87228, 0.24748],
                 [1.00000, 0.83776, 0.20892],
                 [0.99999, 0.80247, 0.17768],
                 [0.99999, 0.76641, 0.15171],
                 [0.99998, 0.72953, 0.12971],
                 [0.99996, 0.69173, 0.11084],
                 [0.99994, 0.65285, 0.09448],
                 [0.99991, 0.61270, 0.08024],
                 [0.99985, 0.57099, 0.06781],
                 [0.99977, 0.52734, 0.05700],
                 [0.99964, 0.48119, 0.04774],
                 [0.99942, 0.43165, 0.04006],
                 [0.99904, 0.37733, 0.03425],
                 [0.99829, 0.31562, 0.03121],
                 [0.99653, 0.24108, 0.03388],
                 [0.98995, 0.14079, 0.05718],
                 [0.95565, 0.04309, 0.16091],
                 [0.89867, 0.01496, 0.28248],
                 [0.83663, 0.00638, 0.38039],
                 [0.77159, 0.00295, 0.45924],
                 [0.70427, 0.00141, 0.52191],
                 [0.63527, 0.00068, 0.56978],
                 [0.56522, 0.00032, 0.60344],
                 [0.49490, 0.00015, 0.62312],
                 [0.42519, 0.00006, 0.62886],
                 [0.35720, 0.00003, 0.62072],
                 [0.29217, 0.00001, 0.59895],
                 [0.23145, 0.00000, 0.56414],
                 [0.17636, 0.00000, 0.51734],
                 [0.12806, 0.00000, 0.46007],
                 [0.08734, 0.00000, 0.39419],
                 [0.05461, 0.00000, 0.32169],
                 [0.02990, 0.00000, 0.24450],
                 [0.01290, 0.00000, 0.16431],
                 [0.00313, 0.00000, 0.08248],
                 [0.00000, 0.00000, 0.00000]]
SUNSET_NSAMPLE = len(SUNSET_COLORS)
SUNSET_NODES = np.linspace(0., 1., SUNSET_NSAMPLE)
SUNSET_CMAP = mcol.LinearSegmentedColormap.from_list("sunset_cmap", list(zip(SUNSET_NODES, SUNSET_COLORS)))
SUNSET_CMAP = SUNSET_CMAP.reversed() ## Reverse the colormap in order to get larger is lighter.
                                     ## Note: Black on white is much more suitable for printing (since it reduces ink bleeding) and more efficient for reading (as reported by many studies) but, unfortunately, when it represents physical quantities it is less easy to interpret because white is generally associated with higher intensities.
SUNSET_CMAP.set_bad('white', 0.) ## Set the color when nan is encountered. Args: (color, opacity).

##cdict = {'red':   [(0.0, 1.0, 1.0),  # red decreases
##                   (1.0, 0.0, 0.0)],
##
##         'green': [(0.0, 0.0, 0.0),  # green increases
##                   (1.0, 1.0, 1.0)],
##
##         'blue':  [(0.0, 0.0, 0.0),  # no blue at all
##                   (1.0, 0.0, 0.0)]}
##
##red_green_cm = LinearSegmentedColormap('RedGreen', cdict, N)


def colormap_to_tikz_code(cmap, nsample, mode):
    """
    Returns a TikZ code version of the given colormap "cmap" using a given number of samples "nsample".
    Example output: 'colormap={temperature}{rgb255=(0,0,128) rgb255=(0,0,255) rgb255=(255,255,255) rgb255=(255,0,0) rgb255=(128,0,0)}'
    """
    ##print(ct.TAG_INFO + "cmap.name =", cmap.name, ", cmap.N =", cmap.N, ", cmap(0.5) =", cmap(0.5))
    
    string = "colormap={" + cmap.name + "}{"
    
    for i in range(nsample):
        val = i/(nsample - 1)
        rgb255 = [int(np.round(255*c)) for c in cmap(val)[0:3]]
        string += "rgb255=(" + str(rgb255[0]) + ","+ str(rgb255[1]) +","+ str(rgb255[2]) +") "
    
    string += "}"
    
    if (mode == "log"):
        string += ", colorbar style={yticklabel={$10^{\\pgfmathprintnumber[fixed relative, precision=3]{\\tick}}$}}"
    
    return string


def array_from_file(field_file, column_name):
    """
    Returns a two-dimensional array containing the function to plot extracted from the column "column_name" in the file "field_file".
    This function also returns the horizontal bounds (xmin, xmax), the vertical bounds (ymin, ymax), and other data extracted from the file.
    """
    try:
        fp = open(field_file, 'r')
    except IOError as e:
        print(ct.TAG_ERROR + "Field file '" + field_file + "' not found, aborting now...")
        return 1
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    data_header = ct.get_header(fp, '%')
    ##holscat = float(ct.get_value_in_string("h/lscat", data_header))  ## Extract the ratio h/lscat.
    fp.close()
    
    ## Find out the bounds of the system:
    point = np.asarray([(int(p['x']), int(p['y'])) for p in data], dtype=int)
    xmin = point[:, 0].min()
    xmax = point[:, 0].max()
    ymin = point[:, 1].min()
    ymax = point[:, 1].max()
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1
    array = np.full((ny, nx), np.nan)
    
    ## Loop on the points to write in 'array':
    for p in data:
        i = ymax - int(p['y'])
        j = int(p['x']) - xmin
        array[i, j] = float(p[column_name])
        ##print(ct.TAG_INFO + "array[", i, ",", j, "] = ", float(p[column_name]))
    
    array = np.ma.array(array, mask=np.isnan(array))  ## Use a mask to escape nan values.
    
    return (array, xmin, xmax, ymin, ymax, data, data_header)


def array_to_tikz(array, xmin, xmax, ymin, ymax, vmin, vmax, unit_length, data, data_header, mode, args, file_path):
    """
    Create a bitmap and a TikZ file importing this bitmap to plot the data.
    This function calls an external script to compile the TikZ file.
    """
    cmap_base = mplt.cm.turbo  ## Use 'turbo' colormap (recommended).
    #cmap = mplt.cm.jet ## Use 'jet' colormap.
    #cmap = SUNSET_CMAP  ## Use custom 'sunset' colormap.
    cmap_color_list = [c for c in cmap_base.colors] ## Extract the colors to resample the colormap.
    cmap_nodes = np.linspace(0., 1., len(cmap_color_list))
    cmap = mcol.LinearSegmentedColormap.from_list("turbo_resampled", list(zip(cmap_nodes, cmap_color_list)), N=1024)
    #mplt.imshow(array, cmap=cmap)  ## Show the plot in live (optional).
    #mplt.colorbar()
    #mplt.show()
    
    norm = mcol.Normalize(vmin=vmin, vmax=vmax)
    image = cmap(norm(array)) ## Create the bitmap image.
    ##print(ct.TAG_INFO + "vmin = ", vmin, ", vmax = ", vmax)
    
    ##pre, ext = os.path.splitext(field_file)
    bitmap_file = file_path + ".png"
    mplt.imsave(bitmap_file, image)  ## Save the raw pixel-constrained bitmap to a PNG file.
    
    tikz_code = """%% Generated on {timestamp} by {my_program} {my_copyright}
{data_header}
\\begin{{tikzpicture}}[%
    /pgfplots/every axis/.append style={{%
        width=0.9\\textwidth,
        {colorbar_string}
    }},
    {boundary_style},
]%
\\begin{{axis}}[%
    title={{{title}}},
    xlabel={{{xlabel}}},
    ylabel={{{ylabel}}},
    xticklabel={{\\pgfmathparse{{{unit_length}*\\tick}}$\\pgfmathprintnumber[fixed relative, precision=3]{{\\pgfmathresult}}$}}, %% Rescale ticks to get x/lscat = (x/h) * (h/lscat), with h/lscat={unit_length}.
    yticklabel={{\\pgfmathparse{{{unit_length}*\\tick}}$\\pgfmathprintnumber[fixed relative, precision=3]{{\\pgfmathresult}}$}},
    colorbar, %% Enable colorbar.
    point meta min={vmin}, %% Set colorbar range.
    point meta max={vmax},
    axis equal image, %% Unit aspect ratio.
    axis line style={{draw=none}}, %% Hide axes.
    tick align=outside,
    enlargelimits={{abs=0.5pt}},
    clip=false, %% Disable axis clipping.
]%
\\addplot graphics[xmin={xmin}, xmax={xmax}, ymin={ymin}, ymax={ymax}]{{{bitmap_file}}};
{boundary_code}
\\end{{axis}}%
\end{{tikzpicture}}%""".format(#
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = ct.MY_COPYRIGHT,
        data_header = data_header,
        title = "\\textbf{Cmd:} \\detokenize{"+ " ".join(args) + "}",
        xlabel = "$x/\\ell$",
        ylabel = "$y/\\ell$",
        unit_length = unit_length,
        colorbar_string = colormap_to_tikz_code(cmap, 41, mode),
        boundary_style = plot_mesh.BOUNDARY_STYLE,
        boundary_code = plot_mesh.boundary_to_tikz_code(data),
        vmin   = vmin,
        vmax   = vmax,
        xmin   = xmin-0.5,
        xmax   = xmax+0.5,
        ymin   = ymin-0.5,
        ymax   = ymax+0.5,
        bitmap_file = "\\jobname.png"
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + ".tikz"
    print(ct.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    ct.compile_tikz(tikz_file) ## Compile the TikZ file.
    return 0


def plot_map(args):
    """
    Plots the given component args[1]='column_name' of the given CSV file args[2]='field_file'.
    The data in the field file are assumed to have the form [x, y, north, south, east, west, f1, f2, ..., fn], where the components 'x' and 'y' are assumed to be integers, and the cardinal directions are the indices of nearest neighbors or boundary conditions.
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 6):
        print(ct.TAG_ERROR + "Invalid number of arguments, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " MODE(lin|log) COLUMN_NAME UNIT_LENGTH VRANGE(auto|vmin:vmax) FIELD_FILE")
        return 1
    
    mode = args[1]  ## Coloring mode ("lin" or "log").
    column_name = args[2]  ## Column in the field file.
    unit_length = float(eval(args[3]))  ## Unit length, typically the best estimate of h/lscat = (L/lscat)/(L/h).
    vrange = args[4]  ## Value range. Either "vmin:vmax" or "auto".
    field_file  = args[5]  ## File containing the field.
    file_path = os.path.splitext(field_file)[0] + "_" + column_name  ## The file path will be used to write new files.
    
    ## Check for possible invalid arguments:
    if (mode != "lin" and mode != "log"):
        print(ct.TAG_ERROR + "Invalid scale mode '" + mode + "', expected 'lin' or 'log', aborting...")
        return 1
    
    if (unit_length <= 0.):
        unit_length = 1.
    
    ## Extract the array to plot from the data file:
    (array, xmin, xmax, ymin, ymax, data, data_header) = array_from_file(field_file, column_name)
    
    if (mode == "log"):
        array = np.log10(array)
    
    ## Find the bounds (vmin, vmax) of the charted function:
    if (vrange == "auto"):
        vmin = array.min() ## If "auto" mode, then extract the depth range of the field [vmin, vmax].
        vmax = array.max()
    elif (":" in vrange):
        vmin, vmax = sorted(list(map(eval, vrange.split(":"))))
    else:
        print(ct.TAG_ERROR + "Invalid value range, expected 'auto' or vmin:vmax (separated by colon), aborting...")
        return 1
    
    ## Plot the array using a bitmap imported in TikZ:
    array_to_tikz(array, xmin, xmax, ymin, ymax, vmin, vmax, unit_length, data, data_header, mode, args, file_path)
    return 0

if (__name__ == '__main__'):
    exit(plot_map(sys.argv))
