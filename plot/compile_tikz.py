#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2024-01-26 at 16:09:23 CET by David Gaspard <david.gaspard@espci.fr>
## Python module providing the compile_tikz() function to compile TikZ files.
## This module is used by many template scripts in the folder ./tikz/
## This file is also callable as a script. When called, it compiles the given TikZ files.
## USAGE : ./plot/compile_tikz.py FILE_1.tikz [ FILE_2.tikz FILE_3.tikz ]
## FILENAME_TIKZ = Path of the TikZ file to be compiled by the present script.
import sys
import os
import shutil

MY_COPYRIGHT = "(c) 2025 David GASPARD (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>"

TAG_INFO = "[INFO] "                     ## Information tag.
TAG_WARN = "[\033[1;93mWARN\033[0m] "    ## Warning tag.
TAG_ERROR = "[\033[1;31mERROR\033[0m] "  ## Error tag.
TAG_EXEC = "[\033[1;95mEXEC\033[0m] "    ## Execution tag.
TAG_USAGE = "[USAGE] "                   ## Usage tag.

LATEX_COMPILER = "pdflatex"    ## Compiler used to compile the LaTeX file.

## Define the LaTeX preamble (do not add LaTeX comments because line breaks are removed):
LATEX_PREAMBLE = r"""\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath, amssymb}
\pagestyle{empty}
\usepackage{pgfplots}
\renewcommand{\Re}{\operatorname{Re}}
\renewcommand{\Im}{\operatorname{Im}}
\newcommand{\D}{\mathop{}\!\mathrm{d}}
\newcommand{\E}{\mathop{}\!\mathrm{e}}
\newcommand{\I}{\mathrm{i}}
\newcommand{\der}[3][]{\frac{\D^{#1} #2}{\D #3^{#1}}}
\newcommand{\pder}[3][]{\frac{\partial^{#1} #2}{\partial #3^{#1}}}
\newcommand{\vect}[1]{\boldsymbol{\mathrm{#1}}}
\newcommand{\matr}[1]{\mathsf{#1}}
\newcommand{\op}[1]{\hat{#1}}
\newcommand{\bra}[1]{\left\langle#1\right|}
\newcommand{\ket}[1]{\left|#1\right\rangle}
\newcommand{\braket}[2]{\left\langle#1\middle|#2\right\rangle}
\newcommand{\avg}[1]{\left\langle#1\right\rangle}
\newcommand{\tavg}[1]{\langle#1\rangle}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\norm}[1]{\left\|#1\right\|}
\newcommand{\tran}[1]{{#1}^{\intercal}} 
\newcommand{\freal}[1]{{#1}_{\rm r}}
\newcommand{\fimag}[1]{{#1}_{\rm i}}
\newcommand{\cc}[1]{{#1}^*}
\newcommand{\herm}[1]{#1^{\dagger}}
\newcommand{\lscat}{\ell_{\rm s}}
\newcommand{\labso}{\ell_{\rm a}}
\definecolor{pyplot1}{HTML}{1F77B4}
\definecolor{pyplot2}{HTML}{FF7F0E}
\definecolor{pyplot3}{HTML}{2CA02C}
\definecolor{pyplot4}{HTML}{D62728}
\usetikzlibrary{decorations.markings}
\tikzset{
	font={\footnotesize},
	ultra thin/.style=  {line width=0.2pt},
	very thin/.style=   {line width=0.4pt},
	thin/.style=        {line width=0.6pt},
	semithick/.style=   {line width=0.8pt},
	thick/.style=       {line width=1.0pt},
	very thick/.style=  {line width=1.4pt},
	ultra thick/.style= {line width=1.8pt},
	every picture/.append style={thin},
	every node/.append style={transform shape},
}
\pgfplotsset{
    compat=1.18,
	every axis/.append style={
        width=0.9\textwidth,
		axis line style={thin},
		every axis title/.style={
            at={(0.5, 1)},
            above,
            align=center,
        },
        legend style={
            thin,
            cells={anchor=west},
        },
        scaled ticks=false,
        ticklabel style={
            /pgf/number format/fixed,
            /pgf/number format/precision=5,
            /pgf/number format/1000 sep={\,},
        },
        yticklabel style={rotate=90},
        enlargelimits=false,
        mark size=1.2,
        axis on top=true,
	},
    every axis plot/.append style={thick},
    table/col sep=comma,
}""".replace('\n', '')

##print(GLOBAL_PREAMBLE_TEX)

def compile_tikz(filename_tikz):
    """
    Compile the TikZ file "filename_tikz" in LaTeX using the global template file "tikz/preamble.tex".
    This function also checks for possible errors. Returns 1 on error, and 0 otherwise. 
    This function assumes that "pdflatex" and "grep" are available.
    """
    ## 1. First check for possible errors:
    if (shutil.which(LATEX_COMPILER) == None):
        print(TAG_ERROR + "LaTeX compiler not found: '" + LATEX_COMPILER + "'...")
        return 1
    if (not os.path.isfile(filename_tikz)):
        print(TAG_ERROR + "TikZ file not found: '" + filename_tikz + "'...")
        return 1
    
    ## 2. Prepare the substitution dictionary:
    jobname = os.path.splitext(filename_tikz)[0]  ## The LaTeX jobname is simply the filename without extension.
    ##search_path = os.path.dirname(filename_tikz)  ## The search path used by PGFPlots to locate data files.
    
    dic = {## Substitution dictionary:
        "compiler": LATEX_COMPILER,
        "preamble": LATEX_PREAMBLE,
        "jobname": jobname,
        "filename_tikz": filename_tikz
    }
    
    ## 3. Prepare the UNIX commands to be executed:
    cmd = "%(compiler)s -jobname '%(jobname)s' '\\documentclass[12pt]{article}%(preamble)s\\begin{document}\\noindent\\input{\detokenize{%(filename_tikz)s}}\\end{document}' | grep -C 1 -wi --color=auto '^!\\|^l\\|error\\|undefined\\|unknown\\|missing\\|runaway\\|misplaced\\|multiply\\|exceeded\\|too\\|ended\\|extra\\|forget\\|forgotten\\|unbounded\\|overfull\\|underfull' " % dic
    
    ##print("COMMAND: " + cmd)
    
    ## 4. Execute the commands:
    print(TAG_INFO + "Compiling TikZ file: '" + filename_tikz + "'...")
    os.system(cmd)
    
    ## 5. Removes LaTeX's auxiliary files:
    flist = [jobname + ".aux", jobname + ".log", jobname + ".out"]
    for f in flist:
        if (os.path.isfile(f)):
            ##print(compile_tikz.TAG_INFO + "Removing file: '" + f + "'...")
            os.remove(f)
    
    return 0

def answer_is_yes(msg):
    """
    Prompts the user with a binary choice.
    Returns True if the answer is "yes", False otherwise.
    """
    while(True):
        ans = input(msg)
        if (ans.lower() in ["y", "yes"]):
            return True
        elif (ans.lower() in ["n", "no"]):
            return False
        print("Please answer yes (y) or no (n)...")

def can_overwrite(filename):
    """
    Prompts the user if the file "filename" can be overwritten.
    Returns True if the file can be overwritten, False otherwise.
    """
    if (os.path.isfile(filename)):
        print(TAG_WARN + "File already exists: '" + filename + "'...")
        if (not answer_is_yes("Overwrite ? (y/n): ")):
            print(TAG_INFO + "OK, keeping file...")
            return False
    return True

def get_header_comment(filename, comment_char):
    """
    Returns the first commented lines of a file assuming the given comment character.
    """
    res = ""
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(comment_char):
                res += line
    return res

def main(args):
    """
    Main function of the script. This function is called when the script is prompted. 
    """
    if (len(args) == 1):
        print(TAG_ERROR + "No input file, doing nothing...")
        print(TAG_USAGE + args[0] + " file_1.tikz [file_2.tikz file_3.tikz ...] ")
        return 1
    
    ## Compile all the given files:
    for f in args[1:]:
        compile_tikz(f)
    
    return 0

if (__name__ == '__main__'):
    exit(main(sys.argv))
