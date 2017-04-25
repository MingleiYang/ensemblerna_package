##########################################################################################################
#Visualization of RNA ensembles
    #getMap
    #getMapDB
    #getRef
    #getRefDB
    #main
##########################################################################################################
import sys
import argparse
import random
from ensemblerna.ErrorCheck import *
from ensemblerna.DBStructs import *
from ensemblerna.DBAnalysis import *
from ensemblerna.PlotVis import *
from ensemblerna.PlotInter import *

##########################################################################################################
#Function to create map of conformational space
#Input: fasta sequence, file header, output directory, map size, ignore flag
#Output: map positions, map 2D matrix, map sequences, map dot bracket structures
##########################################################################################################
def _getMap(fasta, header, dir, size, plotm, rg, rsp, thmax, ignore):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #get map structure
    map = getMapSeqs(dir+'/', fasta, size, rg, rsp, thmax)
    mapseqs = map['mapseqs']
    mapinds = map['mapinds']
    getMapStruct(dir, mapinds, header+'_map.db')
    
    #encode map structure
    bin_mat_2d = encodeStructsNested(dir+header+'_map', rg, ignore)
    map2d = bin_mat_2d['norm2d_']
    mapdb = bin_mat_2d['db_']
    
    #get map frequency
    nestfreq = getNestFreq(map2d)
    mapfreq = nestfreq['freq']
    mapnest = nestfreq['nest']
    mapseqs = nestfreq['seqs']
    maparr = nestfreq['arr']
    
    #plot map
    mappos = plotMap(maparr, mapfreq, mapnest, mapseqs, mapdb, map2d, dir+header+'_map', plotm)
    
    #return map
    return{'mappos':mappos, 'mapnest':mapnest, 'mapseqs':mapseqs, 'mapdb':mapdb, 'map2d':map2d, 'maparr':maparr}

##########################################################################################################
#Function to create map of conformational space from dot-bracket
#Input: dot-bracket file, file header, output directory, map size
#Output: map positions, map 2D matrix, map sequences, map dot bracket structures
##########################################################################################################
def _getMapDB(mapdb, header, dir, size, plotm, rg, ignore):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #get map structure
    cmd = 'cat ' + mapdb + ' > ' + dir+header+'_map.db'
    subprocess.check_output(cmd, shell=True)
    
    #encode map structure
    bin_mat_2d = encodeStructsNested(dir+header+'_map', rg, ignore)
    map2d = bin_mat_2d['norm2d_']
    mapdb = bin_mat_2d['db_']
    
    #get map frequency
    nestfreq = getNestFreq(map2d)
    mapfreq = nestfreq['freq']
    mapnest = nestfreq['nest']
    mapseqs = nestfreq['seqs']
    maparr = nestfreq['arr']
    
    #plot map
    mappos = plotMap(maparr, mapfreq, mapnest, mapseqs, mapdb, map2d, dir+header+'_map', plotm)
    
    #return map
    return{'mappos':mappos, 'mapnest':mapnest, 'mapseqs':mapseqs, 'mapdb':mapdb, 'map2d':map2d, 'maparr':maparr}


##########################################################################################################
#Function to create reference visualization
#Input: map object, fasta sequence, file header,output directory, ignore flag, number of samples, shape data
#Output: csv, db, pdf, png
##########################################################################################################
def _getRef(map, fasta, header, dir, rg, plotint, ignore, numsamp, shape=None):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #get structure
    if shape is None:
        getStruct(dir, header+'.db', numsamp)
    else:
        getStructSHAPE(dir, header+'.db', numsamp)

    #read structure
    bin_mat_2d = encodeStructsNested(dir+header, rg, ignore)
    ref2d = bin_mat_2d['norm2d_']
    refdb = bin_mat_2d['db_']
    
    #get reference frequency
    nestfreq = getNestFreq(ref2d)
    reffreq = nestfreq['freq']
    refnest = nestfreq['nest']
    refseqs = nestfreq['seqs']
    refarr = nestfreq['arr']

    #get map
    mappos = map['mappos']
    mapnest = map['mapnest']
    mapseqs = map['mapseqs']
    mapdb = map['mapdb']
    map2d = map['map2d']
    maparr = map['maparr']
    
    #plot reference
    ref = plotRef(reffreq, refnest, refarr, mappos, mapnest, mapseqs, mapdb, maparr, map2d, dir+header, refdb, refseqs, rg)

    #plot interactive
    if plotint == 'T':
        structs = ref['structs']
        diversity = ref['diversity']
        freq = ref['freq']
        plotInteractive(dir, header, mappos, freq, structs, diversity, maparr, rg)

    return True

##########################################################################################################
#Function to create reference visualization with dot bracket
#Input: map object, dot bracket, file header, ignore flag, output directory
#Output: csv, db, pdf, png
##########################################################################################################
def _getRefDB(map, db, header, dir, rg, plotint, ignore, md=None):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #remove map temporary files
    if md is None:
        cmd = 'rm ' + dir + '/temp*'
        subprocess.check_output(cmd, shell=True)
    
    #get reference structure
    cmd = 'cat ' + db + ' > ' + dir+header+'.db'
    subprocess.check_output(cmd, shell=True)

    #read structure
    bin_mat_2d = encodeStructsNested(dir+header, rg, ignore)
    ref2d = bin_mat_2d['norm2d_']
    refdb = bin_mat_2d['db_']
    
    #get reference frequency
    nestfreq = getNestFreq(ref2d)
    reffreq = nestfreq['freq']
    refnest = nestfreq['nest']
    refseqs = nestfreq['seqs']
    refarr = nestfreq['arr']
    
    #get map
    mappos = map['mappos']
    mapnest = map['mapnest']
    mapseqs = map['mapseqs']
    mapdb = map['mapdb']
    map2d = map['map2d']
    maparr = map['maparr']
    
    #plot reference
    ref = plotRef(reffreq, refnest, refarr, mappos, mapnest, mapseqs, mapdb, maparr, map2d, dir+header, refdb, refseqs, rg)

    #plot interactive
    if plotint == 'T':
        structs = ref['structs']
        diversity = ref['diversity']
        freq = ref['freq']
        plotInteractive(dir, header, mappos, freq, structs, diversity, maparr, rg)
    
    return True


#runs main part of script
def main():
    
    #check programs loaded
    checkRNAStructCMD()
    
    #set seed
    random.seed(113)
    
    #parse command line
    parser = argparse.ArgumentParser(prog='EnsembleRNA', usage='ensemblerna <fasta file> <output directory> [options]', description='Visualize the structural ensemble for a given RNA. For more information please see the README file, the Documentation file, or visit http://ribosnitch-ensemblerna.rhcloud.com')
    parser.add_argument('-v', '--version',  action='version', version='%(prog)s 1.0')
    parser.add_argument('input', metavar='fasta', type=str, help='Reference fasta file. Maximum sequence length is 2500 nucleotides (Required)')
    parser.add_argument('output', metavar='outdir', type=str, help='Output directory (Required)')
    parser.add_argument('-sh', metavar='--shape', type=str, dest='shape_file', default=None, help='Includes shape data in the reference ensemble prediction. Ignored if -d flag is used (Default is None)')
    parser.add_argument('-d', metavar='--db', type=str, dest='db_file', default=None, help='Dot-bracket structures for reference ensemble (Default is None)')
    parser.add_argument('-m', metavar='--map', type=str, dest='map_file', default=None, help='Sequence to create the map of conformational space. Ignored if -md flag is used (Default is reference fasta file)')
    parser.add_argument('-md', metavar='--mapdb', type=str, dest='map_dbfile', default=None, help='Dot-bracket structures for the map of conformational space. A previously created map can be used to project new ensembles onto the same space. (Default is None)')
    parser.add_argument('-s', metavar='--size', dest='size', type=int, default=10, help='Number of sequences for the map of conformational space. Higher numbers increase structural diversity. Ignored if -md flag is used (Default is 10)')
    parser.add_argument('-p', metavar='--plotmap', dest='plotm', type=str, choices=['T', 'F'], default='T', help='Plot the map T/F (Default is T)')
    parser.add_argument('-r', metavar='--range', dest='nucrg', type=int, nargs=2, default=None, help='Range of nucleotides to visualize. Predicted structures will still include the full length of the input RNA, but only the given range will be plotted (Default is 1 to sequence length)')
    parser.add_argument('-maxd', metavar='--maxdistance', dest='maxd', type=int, default=None, help='Maximum number of bases between the two nucleotides in a pair (Default is no restriction)')
    parser.add_argument('-t', metavar='--temperature', dest='tempcalc', type=int, default=None, help='Temperature at which the calculation takes place in Kelvin (Default is 310.15 K)')
    parser.add_argument('-si', metavar='--SHAPEintercept', dest='sint', type=int, default=None, help='Intercept used with SHAPE restraints. Ignored if -d flag is used (Default is -0.6 kcal/mol)')
    parser.add_argument('-sm', metavar='--SHAPEslope', dest='slope', type=int, default=None, help='Slope used with SHAPE restraints. Ignored if -d flag is used (Default is 1.8 kcal/mol)')
    parser.add_argument('-pi', metavar='--plotinteractive', dest='plotint', type=str, choices=['T', 'F'], default='T', help='Plot the interactive file T/F (Default is T)')
    parser.add_argument('-th', metavar='--threadmax', dest='thmax', type=int, default=1, help='Maximum number of threads for multi-threading. (Default is 1)')
    parser.add_argument('-i', metavar='--ignorestems', dest='ignore', type=int, default=3, help='Ignore stems with fewer than i base pairs. (Default is 3)')
    parser.add_argument('-n', metavar='--num', dest='numsamp', type=int, default=1000, help='Number of Boltzmann sampled structures to produce for the visualization. (Default is 1000)')
    args = parser.parse_args()
    fasta_file = args.input
    outdir = args.output
    shape_file = args.shape_file
    map_file = args.map_file
    size = args.size
    db_file = args.db_file
    map_dbfile = args.map_dbfile
    plotm = args.plotm
    nucrg = args.nucrg
    maxd = args.maxd
    sint = args.sint
    slope = args.slope
    tempcalc = args.tempcalc
    plotint = args.plotint
    thmax = args.thmax
    ignore = args.ignore
    numsamp = args.numsamp
    
    #initialize variables
    header = fasta_file.split('/')[-1]
    header = header.split('.')[0]

    #edit outdir
    if(outdir[-1] == '/'):
        outdir = outdir[:-1]

    #read fasta file
    fasta = checkFasta(fasta_file)

    #check inputs for errors
    outdir = checkDir(outdir)
    if nucrg == None:
        rg = range(0,len(fasta),1)
    else:
        rg = checkRange(nucrg, len(fasta))
    checkSize(size, len(fasta), rg)
    checkIgnore(ignore, len(fasta), rg)
    checkThmax(thmax)
    checkNumsamp(numsamp)

    #check RNAstructure parameters
    rsp = ''
    ssp = ''
    if maxd is not None:
        rsp = rsp+' -md '+str(maxd)
    if tempcalc is not None:
        rsp = rsp+' -t '+str(tempcalc)
    if sint is not None:
        ssp = ssp+' -si '+str(sint)
    if sint is not None:
        ssp = ssp+' -sm '+str(slope)

    #check reference
    if db_file is None:
        #choose shape function (shape or none)
        if shape_file is None:
            shape = None
            checkRNAStruct(fasta, outdir, rsp)
        else:
            shape = 1
            checkSHAPE(shape_file, outdir)
            checkRNAStructSHAPE(fasta, outdir, rsp, ssp)
    else:
        #check dot bracket
        checkDB(db_file, outdir, len(fasta))

    #choose map function (db or fasta)
    if map_dbfile is not None:
        #check dot bracket
        checkDB(map_dbfile, outdir, len(fasta))
        
        #get map
        map = _getMapDB(map_dbfile, header, outdir, size, plotm, rg)

    else:
        #read map file
        if map_file is not None:
            mfasta = checkFasta(map_file, len(fasta))
        else:
            mfasta = fasta

        #get map
        map = _getMap(mfasta, header, outdir, size, plotm, rg, rsp, thmax, ignore)

    #choose reference function (db or fasta)
    if db_file is not None:
        #visualize reference
        if map_dbfile is None:
            ref = _getRefDB(map, db_file, header, outdir, rg, plotint, ignore, numsamp)
        else:
            ref = _getRefDB(map, db_file, header, outdir, rg, plotint, ignore, md=1)
    else:
        #visualize reference
        ref = _getRef(map, fasta, header, outdir, rg, plotint, ignore, numsamp, shape)


