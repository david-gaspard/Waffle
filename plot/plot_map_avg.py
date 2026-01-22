#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2026-01-22 at 12:44:51 CET by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot a scalar field defined over a square lattice with holes. The input data must have the form [x, y, f1(x,y), f2(x,y), ..., fn(x,y)].
import sys, os, datetime, csv
import numpy as np
import compile_tikz as ct
import plot_map


def avg_array_from_file(field_file, ncol):
    """
    Returns a two-dimensional array containing the function to plot computed from the average of the first "ncol" columns in the file "field_file".
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
    
    ## Check for possible invalid number of columns:
    ncolmax = len(data[0])-6
    if (ncol <= 0 or ncol > ncolmax):
        print(ct.TAG_ERROR + "Invalid number of columns " + str(ncol) + ", expected in 1.." + str(ncolmax) + ", aborting...")
        return 1
    
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
        array[i, j] = 0.
        for k in range(ncol):
            colname = "I" + str(k)
            array[i, j] += float(p[colname])
        array[i, j] /= ncol
    
    array = np.ma.array(array, mask=np.isnan(array))  ## Use a mask to escape nan values.
    
    return (array, xmin, xmax, ymin, ymax, data, data_header)


def plot_map_avg(args):
    """
    Plots the average of the first "ncol" functions of the given CSV file args[2]='field_file'.
    The data in the field file are assumed to have the form [x, y, north, south, east, west, f1, f2, ..., fn], where the components 'x' and 'y' are assumed to be integers, and the cardinal directions are the indices of nearest neighbors or boundary conditions.
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 6):
        print(ct.TAG_ERROR + "Invalid number of arguments, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " MODE(lin|log) NCOL UNIT_LENGTH VRANGE(auto|vmin:vmax) FIELD_FILE")
        return 1
    
    mode = args[1]  ## Coloring mode ("lin" or "log").
    ncol = int(args[2])  ## Number of columns to average.
    unit_length = float(eval(args[3]))  ## Unit length, typically the best estimate of h/lscat = (L/lscat)/(L/h).
    vrange = args[4]  ## Value range. Either "vmin:vmax" or "auto".
    field_file  = args[5]  ## File containing the field.
    file_path = os.path.splitext(field_file)[0] + "_avg" + str(ncol)  ## The file path will be used to write new files.
    
    ## Check for possible invalid arguments:
    if (mode != "lin" and mode != "log"):
        print(ct.TAG_ERROR + "Invalid scale mode '" + mode + "', expected 'lin' or 'log', aborting...")
        return 1
    
    if (unit_length <= 0.):
        unit_length = 1.
    
    ## Extract the array to plot from the data file:
    (array, xmin, xmax, ymin, ymax, data, data_header) = avg_array_from_file(field_file, ncol)
    
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
    plot_map.array_to_tikz(array, xmin, xmax, ymin, ymax, vmin, vmax, unit_length, data, data_header, mode, args, file_path)
    return 0


if (__name__ == '__main__'):
    exit(plot_map_avg(sys.argv))
