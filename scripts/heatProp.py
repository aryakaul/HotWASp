import argparse
import sys
import numpy as np

def parseargs(arguments):
    parser = argparse.ArgumentParser(description="Check help flag")
    parser.add_argument("-n", "--network", type=argparse.FileType('r', encoding='UTF-8'), required=True, help= "Direct path to a gene,gene unweighted network edge list. File format specified in README")
    parser.add_argument("-t", "--uppertimebound", required=False, help="Upper limit to the maximum time of t", default=3, type=int)
    parser.add_argument("-g", "--gwas", required=True, help="Direct path to the gwas results being boosted. File format specified in README", type=argparse.FileType('r', encoding='UTF-8'))
    parser.add_argument("-o", "--output", required=False, help="Full path to output file. Default is current directory in file 'HotWASp-results.txt'", default="./HotWASp-results.txt")
    return parser.parse_args()

def readNetwork(networkPath):
    with open(networkPath, 'r') as filein:
        for lines in filein:

def readGWAS(gwasPath):
    with open(gwasPath, 'r') as filein:
        for lines in filein:



def heatProp (adjMatrix, heat, prop = 2, time = 1):
    degree = np.sum(adjMatrix,1)
    heatProp = np.identity(adjMatrix.shape[0])
    laplace = heatProp
    print("Calculating Laplace")
    for i in range(len(degree)):
        laplace[i,i] = degree[i]
    laplace -= adjMatrix
    laplace = laplace/time
    lpTest = np.sum(laplace,1)
    lpNeg = []
    for i in lpTest:
        if i < 0:
            lpNeg.append(i)
    print (lpNeg)
    laplacePower = laplace
    p = 1
    while prop > 0:
        print("Calculating Propagation %i"%p)
        heatProp += laplacePower/ math.factorial(p)
        p += 1
        prop -= 1
        if prop != 0:
            laplacePower = np.matmul(laplacePower,laplace)
    print("Calculating new heat")
    newHeat = np.matmul(heat,heatProp)
    return newHeat

def output(heatRanks, output):
    with open(output, 'w') as fileout:
        for tups in heatRanks:
            gene = tup[0]
            heat = tup[1]
            fileout.write("%s\t%s" % (gene, heat))

def maxheats(masterheatRanks, output):
    return finalheats

def main():
    args = parseargs(sys.argv[1:])
    adjMatrix = readNetwork(args.network)
    heat = readGWAS(args.gwas)
    timeVals = np.linspace(0, args.uppertimebound, 10)
    masterheatranks = []
    for t in timeVals:
        heatranks = heatProp(adjMatrix, heat, t)
        masterheatranks.append(heatranks)
    finalheats = maxheats(masterheatranks)
    output(finalheats, args.output)

if __name__=="__main__":
	main()
