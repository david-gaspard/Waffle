#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2026-01-22 at 12:44:51 CET by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot a scalar field defined over a square lattice with holes. The input data must have the form [x, y, f1(x,y), f2(x,y), ..., fn(x,y)].
import sys, os, argparse, datetime, csv
import numpy as np
import compile_tikz as ct
import plot_map


def avg_array_from_file(field_file, nfield):
    """
    Returns a two-dimensional array containing the function to plot computed from the average of the first "nfield" columns in the file "field_file".
    This function also returns the horizontal bounds (xmin, xmax), the vertical bounds (ymin, ymax), and other data extracted from the file.
    """
    try:
        fp = open(field_file, 'r')
    except IOError as exc:
        raise IOError("Failed to open field file '" + field_file + "'.") from exc
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    data_header = ct.get_header(fp, '%')
    ##holscat = float(ct.get_value_in_string("h/lscat", data_header))  ## Extract the ratio h/lscat.
    fp.close()
    
    ## Check for possible invalid number of columns:
    nfieldmax = len(data[0])-6
    if (nfield <= 0 or nfield > nfieldmax):
        print(ct.TAG_ERROR + "Invalid number of columns " + str(nfield) + ", expected in 1.." + str(nfieldmax) + ", aborting...")
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
    
    ## Find out the field names (detecting those with numbers):
    fieldnames = [k for k, v in data[0].items() if k[-1].isdigit()]
    
    ## Loop on the points to write in 'array':
    for p in data:
        i = ymax - int(p['y'])
        j = int(p['x']) - xmin
        array[i, j] = 0.
        for k in range(nfield):
            colname = fieldnames[k]
            array[i, j] += float(p[colname])
        array[i, j] /= nfield
    
    array = np.ma.array(array, mask=np.isnan(array))  ## Use a mask to escape nan values.
    
    return (array, xmin, xmax, ymin, ymax, data, data_header)


def plot_map_avg(args):
    """
    Plots the average of the first "nfield" functions of the given CSV file args[2]='field_file'.
    The data in the field file are assumed to have the form [x, y, north, south, east, west, f1, f2, ..., fn], where the components 'x' and 'y' are assumed to be integers, and the cardinal directions are the indices of nearest neighbors or boundary conditions.
    """
    ## Check for possible invalid arguments:
    available_modes = ["lin", "log", "asinh", "sqrt"]
    if (args.mode not in available_modes):
        print(ct.TAG_ERROR + "Invalid scale mode '" + args.mode + "', expected in " + str(available_modes) + ", aborting...")
        return 1
    
    if (args.unit <= 0.):
        print(ct.TAG_ERROR + "Unit length cannot be negative or zero (received " + args.unit + "), aborting...")
        return 1
    
    ## Extract the array to plot from the data file:
    try:
        (array, xmin, xmax, ymin, ymax, data, data_header) = avg_array_from_file(args.field_file, args.n_field)
    except IOError as exc:
        print(ct.TAG_ERROR + "Field file '" + args.field_file + "' not found, aborting now...")
        return 1
    
    if (args.mode == "log"):
        array = np.log10(array)
    elif (args.mode == "asinh"):
        array = np.arcsinh(array)
    elif (args.mode == "sqrt"):
        array = np.sqrt(array)
    
    ## Find the bounds (vmin, vmax) of the charted function:
    if (args.vrange == "auto"):
        vmin = array.min() ## If "auto" mode, then extract the depth range of the field [vmin, vmax].
        vmax = array.max()
    elif (":" in args.vrange):
        vmin, vmax = sorted(list(map(float, args.vrange.split(":"))))
    else:
        print(ct.TAG_ERROR + "Invalid value range, expected 'auto' or 'vmin:vmax' (separated by colon), aborting...")
        return 1
    
    ## Plot the array using a bitmap imported in TikZ:
    file_path = os.path.splitext(args.field_file)[0] + "_avg" + str(args.n_field)  ## The file path will be used to write new files.
    plot_map.array_to_tikz(array, xmin, xmax, ymin, ymax, vmin, vmax, data, data_header, args, file_path)
    return 0

def parseArguments():
    """
    Parse the arguments of the program.
    """
    parser = argparse.ArgumentParser()  ## Create argument parser.
    
    ## Positional mandatory arguments:
    parser.add_argument("n_field", type=int, help="Number of fields over which to average.")
    parser.add_argument("field_file",  type=str, help="Full file path of the field file.")
    
    ## Optional arguments:
    parser.add_argument("-m", "--mode", type=str, default="lin",
        help="Scale mode, either 'lin', 'log', 'asinh', or 'sqrt'.")
    parser.add_argument("-u", "--unit", type=float, default=1.,
        help="Unit length. Two common choices are either the value of 'h/lscat' or 'kh'.")
    parser.add_argument("-v", "--vrange", type=str, default="auto",
        help="Range of values of the map. Either 'auto' or 'vmin:vmax'.")
    parser.add_argument("-e", "--epsilon", type=float, default=1.3,
        help="Tolerance of the Ramer-Douglas-Peucker algorithm. Larger is simpler (eps=0 to disable).")
    
    return parser.parse_args() ## Return the parsed arguments.


if (__name__ == '__main__'):
    args = parseArguments()  ## Parse the arguments
    exit(plot_map_avg(args))
