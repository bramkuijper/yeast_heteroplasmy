#!/usr/bin/env python

import os, re, sys

first = True

# process all the data and sort out what we want
def process_data(linelist):

    meanp = [0.0,0.0,0.0,0.0]
    ssp = [0.0,0.0,0.0,0.0]

    for line in linelist:
        vals = line.strip().split(";")
        
        if int(vals[0]) != 25:
            continue

        spore = int(vals[2])
        freq_c1 = float(vals[3])
        freq_pop = float(vals[4]);

        meanp[spore] += freq_c1 * freq_pop
        ssp[spore] += freq_c1 * freq_c1 * freq_pop
    

    for i in range(0,4):
        ssp[i] -= meanp[i]*meanp[i]

    return(meanp + ssp)




def analyze_parameters(lines,first=False):

    pars = {}

    for line in lines:
        mobj = re.match("(.*);(.*)",line)
        if mobj != None:
            pars[mobj.group(1)] = mobj.group(2)

    return(pars)

def analyze_file(filename):

    global first;

    # open file; read data
    f = open(filename)
    fl = f.readlines()
    f.close

    if len(fl) < 3:
        return

    flhead = fl[0]

    lc = len(fl)
    parline = -1

    linerange = range(0,lc)
    linerange.reverse()

    # search the parameter line
    for lineno in linerange:

        if re.match("^type",fl[lineno]) != None:
            parline = lineno
            break

    if parline == -1:
        return

    parameters = analyze_parameters(fl[parline:])

    datvals = process_data(fl[1:parline-2])

    flhead = "spore1;spore2;spore3;spore4;var_spore1;var_spore2;var_spore3;var_spore4;"

    if first:
        print ";".join(parameters.keys()) + ";" + flhead + "file"
        first = False

    print ";".join(parameters.values()) + ";" + ";".join([ str(i) for i in datvals]) + ";" + filename


def visit(arg, dirname, names):
    for name in names:
        if re.match("sim_*",name) != None:
            data = analyze_file(dirname + "/" + name)



os.path.walk(sys.argv[1], visit, None)
