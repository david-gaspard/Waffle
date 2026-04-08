#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-07-30 at 11:55:37 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Based on a previous script created on 2025-07-20 at 19:15:01 CEST by the same author.
## Python script to plot a mesh file. This version optimizes the TikZ code to create paths of one piece.
## This optimization allows to use the path decoration in TikZ, hence improving the rendering.
import sys, os, datetime, csv
import numpy as np
import compile_tikz as ct

BOUNDARY_STYLE_OLD = """mirror/.style={black},
    input/.style={red, decoration={markings, mark=between positions 0 and 1 step 1.3mm with {\\fill (-0.5mm, 1mm) -- (0mm, 0mm) -- (0.5mm, 1mm) --cycle;}, pre length=0.5mm, post length=0.5mm}, postaction={decorate}},
    output/.style={blue, decoration={markings, mark=between positions 0 and 1 step 1.3mm with {\\fill (-0.5mm, 0mm) -- (0mm, 1mm) -- (0.5mm, 0mm) --cycle;}, pre length=0.5mm, post length=0.5mm}, postaction={decorate}},
    open/.style={green!80!black, decoration={markings, mark=between positions 0 and 1 step 1.3mm with {\\fill (-0.5mm, 0mm) arc[start angle=180, end angle=0, radius=0.5mm] --cycle;}, pre length=0.5mm, post length=0.5mm}, postaction={decorate}}"""

BOUNDARY_STYLE = """input/.style={%% #1=Spacing between arrows.
        red,
        thick,
        decoration={%
            transform={scale=0.8},
            markings,
            mark=between positions 0 and 1 step #1 with {%
                \\fill (-0.7mm, 1.6mm) -- (0mm, 0.2mm) -- (0.7mm, 1.6mm) -| (0.2mm, 2.6mm) -| (-0.2mm, 1.6mm) --cycle;
            },
            pre length={#1/2},
        },
        postaction={decorate}
    },
    input/.default={3mm}, %% If no argument is given then this kicks in as argument.
    output/.style={%% #1=Spacing between arrows.
        blue,
        thick,
        decoration={%
            transform={scale=0.8},
            markings,
            mark=between positions 0 and 1 step #1 with {%
                \\fill (-0.2mm, 0mm) |- (-0.7mm, 1mm) -- (0mm, 2.4mm) -- (0.7mm, 1mm) -| (0.2mm, 0mm) --cycle;
            },
            pre length={#1/2},
        },
        postaction={decorate}
    },
    output/.default={3mm}, %% If no argument is given then this kicks in as argument.
    mirror/.style={%
        black,
        thick,
        decoration={%
            transform={scale=0.8},
            markings,
            mark=between positions 0 and 1 step 0.7mm with {%
                \\draw[thin] (0mm, 0mm) -- (0.5mm, 1mm);
            },
        },
        postaction={decorate}
    },
    open/.style={draw=none}"""


def point_equal(p1, p2):
    """
    Returns true if the two points are equal up to a small number.
    """
    return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1]) < 1e-8

def merge_segments(segments):
    """
    Merge the segments using a single-pass algorithm and return the number of merges.
    This functions must be called iteratively until the number of merges is zero (i.e., the list does not change anymore).
    """
    copied = len(segments) * [False]  ## Mask indicating if the segment has been treated.
    nmerge = 0  ## Total number of merges.
    
    ##print(ct.TAG_INFO + "Number of segments before merge =", len(segments))
    
    ## 1. Merge segments:
    for i in range(len(segments)):
        if (not copied[i]):
            for j in range(len(segments)):
                if (j != i and not copied[j] and segments[i][0] == segments[j][0] and point_equal(segments[i][1][-1], segments[j][1][0])): 
                    ## If same kind of boundary condition and end points matches, then merge.
                    segments[i][1] = np.concatenate((segments[i][1], segments[j][1][1:]))
                    copied[j] = True
                    nmerge += 1
    
    ## 2. Remove copied segments:
    for i in range(len(segments)-1, -1, -1):
        if (copied[i]):
            del segments[i]
    
    ##print(ct.TAG_INFO + "Number of merges =", nmerge)
    
    return nmerge

def distance_ortho(a, b, p):
    """
    Compute the orthogonal distance of point "p" with respect to the line segment (a, b).
    This is the standard distance used in the Douglas-Peucker algorithm.
    Note that "a", "b", and "p" are assumed to be numpy arrays (2D or 3D).
    """
    u = 0.
    dab2 = np.dot(b-a, b-a)
    if (dab2 > 1.0e-30):  ## Avoid dividing by zero when points are the same (a=b).
        u = np.dot(p-a, b-a)/dab2
    
    if (u <= 0.):
        d = np.linalg.norm(p - a)
    elif (u >= 1.):
        d = np.linalg.norm(p - b)
    else:
        d = np.linalg.norm(p - (a + (b-a)*u))
    
    return d

def find_most_distant(points, indices, epsilon):
    """
    Find the most distant points from the subpath of "points" defined by the list of indices "indices".
    """
    nadded = 0  ## Number of added points to the list "indices".
    
    for iseg in range(len(indices)-1):  ## Loop over segments of the subpath.
        ## Find the most distant point from the segment:
        dmax = 0.
        imax = 0.
        for ip in range(indices[iseg]+1, indices[iseg+1]):
            d = distance_ortho(points[indices[iseg]], points[indices[iseg+1]], points[ip])
            if (d > dmax):
                dmax = d
                imax = ip
        
        ## If the most distant point is farther than epsilon, then save it:
        if (dmax > epsilon):
            indices.append(imax)
            nadded += 1
    
    indices.sort()  ## Sort the list of saved indices.
    ##print("[INFO] Nadded = ", nadded)
    return nadded != 0

def simplify_path_rdp(points, epsilon):
    """
    Simplify the path defined by "points" (list of 2D or 3D numpy arrays) using the Ramer-Douglas-Peucker algorithm.
    The tolerance "epsilon" must be interpreted as the maximum suppression distance for points measured in terms of the function "dist(a, b, p)".
    Greater "epsilon" makes more simplified paths.
    """
    n = points.shape[0]  ## Get the total number of points.
    indices = [0, n-1]  ## List of retained points.
    
    while True:
        if not find_most_distant(points, indices, epsilon):
            break
    
    return points[indices]


def boundary_to_tikz_code(data):
    """
    Returns a TikZ code version of the given 'data', a list of dictionaries generated by csv.DictReader().
    The data contains the fields [x, y, north, south, east, west]. The components 'x' and 'y' are assumed to be integers, 
    and the direction components can be either point indices or boundary conditions ('mirror', 'open', 'input', or 'output').
    """
    ## 1. Create a list of segments on the boundary. The boundary must be travelled in the clockwise direction:
    segments = []  ## Format: [["mirror", ((x1, y1), (x2, y2))], ["open", ((x1, y1), (x2, y2))], ...]
    for p in data:
        if (not p['north'].isdigit()):
            x1 = int(p['x']) - 0.5
            y1 = int(p['y']) + 0.5
            x2 = int(p['x']) + 0.5
            y2 = int(p['y']) + 0.5
            segments.append([p['north'], np.array([(x1, y1), (x2, y2)])])
        if (not p['south'].isdigit()):
            x1 = int(p['x']) + 0.5
            y1 = int(p['y']) - 0.5
            x2 = int(p['x']) - 0.5
            y2 = int(p['y']) - 0.5
            segments.append([p['south'], np.array([(x1, y1), (x2, y2)])])
        if (not p['east'].isdigit()):
            x1 = int(p['x']) + 0.5
            y1 = int(p['y']) + 0.5
            x2 = int(p['x']) + 0.5
            y2 = int(p['y']) - 0.5
            segments.append([p['east'], np.array([(x1, y1), (x2, y2)])])
        if (not p['west'].isdigit()):
            x1 = int(p['x']) - 0.5
            y1 = int(p['y']) - 0.5
            x2 = int(p['x']) - 0.5
            y2 = int(p['y']) + 0.5
            segments.append([p['west'], np.array([(x1, y1), (x2, y2)])])
    
    ## 2. Merge the segments iteratively:
    while True:
        nmerge = merge_segments(segments)
        if (nmerge == 0):
            break
    
    ## 3. Simplify the path:
    epsilon = 1.3  ## Tolerance of the Ramer-Douglas-Peucker algorithm (value eps=1.3 is ok).
    for seg in segments:
        seg_simple = simplify_path_rdp(seg[1], epsilon)
        ##print("[INFO] Simplification ", seg[1].shape[0], " -> ", seg_simple.shape[0])
        seg[1] = seg_simple
    
    ## 4. Convert the path to TikZ code:
    string = "\\begin{scope}%% Draw boundaries\n"
    
    for seg in segments:
        npoint = seg[1].shape[0]
        closed = point_equal(seg[1][0], seg[1][-1])
        if (closed):   ## If the path is closed, then ignore the last point.
            npoint -= 1
        string += "\\draw[{bnd}] (axis cs:{x}, {y})".format(bnd=seg[0], x=seg[1][0,0], y=seg[1][0,1])
        for i in range(1, npoint): ## Loop on the points, excluding the first one.
            string += " -- (axis cs:{x}, {y})".format(x=seg[1][i,0], y=seg[1][i,1])
            if (i%6 == 5 and i != npoint-1): ## Avoid too long lines.
                string += "\n"
        if (closed):
            string += " -- cycle"
        string += ";\n"
    
    return string + "\\end{scope}%"

def plot_mesh(args):
    """
    Generate Tikz code to plot the mesh file (CSV format) given in the first argument: args[1]='filename'.
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 2):
        print(ct.TAG_ERROR + "No input file, doing nothing...")
        print(ct.TAG_USAGE + args[0] + " MESH_FILE")
        return 1
    
    mesh_file = args[1]
    file_path = os.path.splitext(mesh_file)[0]  ## The file path is the filename without its extension (used to write new files). 
    
    try:
        fp = open(mesh_file, 'r')
    except IOError as e:
        print(ct.TAG_ERROR + "Field file '" + mesh_file + " not found, aborting now...")
        return 1
    
    data = list(csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True))
    data_header = ct.get_header(fp, '%')
    
    ## Extract the bounds of the mesh:
    xmesh = np.asarray([int(p['x']) for p in data], dtype=int)
    ymesh = np.asarray([int(p['y']) for p in data], dtype=int)
    xmin = xmesh.min()
    xmax = xmesh.max()
    ymin = ymesh.min()
    ymax = ymesh.max()
    
    tikz_code = """%% Generated on {timestamp} by {my_program} {my_copyright}
{data_header}
\\begin{{tikzpicture}}[%
    {boundary_style},
]%
\\begin{{axis}}[%
    title={{{title}}},
    xlabel={{{xlabel}}},
    ylabel={{{ylabel}}},
    xmin={xmin}, xmax={xmax},
    ymin={ymin}, ymax={ymax},
    axis equal,
    enlargelimits=true, %% Allow for larger view.
]%
{boundary_code}
\\end{{axis}}%
\\end{{tikzpicture}}%""".format(
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = ct.MY_COPYRIGHT,
        data_header = data_header,
        title = "\\textbf{Cmd:} \\detokenize{"+ " ".join(args) + "}",
        xlabel = "$x$",
        ylabel = "$y$",
        xmin   = xmin-0.5,
        xmax   = xmax+0.5,
        ymin   = ymin-0.5,
        ymax   = ymax+0.5,
        boundary_style = BOUNDARY_STYLE,
        boundary_code = boundary_to_tikz_code(data)
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + '.tikz'
    print(ct.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    ct.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    fp.close()
    return 0

if (__name__ == '__main__'):
    exit(plot_mesh(sys.argv))
