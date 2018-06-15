import argparse
import sys
import numpy as np
import math

def parseargs(arguments):
    parser = argparse.ArgumentParser(description="Check help flag")
    parser.add_argument("-n", "--network", type=argparse.FileType('r', encoding='UTF-8'), required=True, help= "Direct path to a gene,gene unweighted network edge list. File format specified in README")
    parser.add_argument("-t", "--uppertimebound", required=False, help="Upper limit to the maximum time of t", default=3, type=int)
    parser.add_argument("-g", "--gwas", required=True, help="Direct path to the gwas results being boosted. File format specified in README", type=argparse.FileType('r', encoding='UTF-8'))
    parser.add_argument("-o", "--output", required=False, help="Full path to output file. Default is current directory in file 'HotWASp-results.txt'", default="./HotWASp-results.txt")
    return parser.parse_args()

def readNetwork(networkPath, header = True,sep = "\t"):
    netList = []
    with open(networkPath, 'r') as filein:
        if header == True:
            netList = [x.strip().split(sep) for x in filein.readlines()[1:]]
        else: 
            netList = [x.strip().split(sep) for x in filein.readlines()]
        
        
    source = [x[0] for x in netList]
    target = [x[1] for x in netList]
    
    geneList = list(set(source) | set(target))
    geneList.sort()
    
    del source
    del target
    
    length = len(netList)
    count = 0
    netDict = {}
    
    print("Dict Making")
    while len(netList) > 0:
        i = netList.pop()
        if count% 5000 == 0:
            print ("...%i/%i"%(count,length))
        count +=1
        if i[0] not in netDict:
            netDict[i[0]] = set()
            
        if i[1] not in netDict:
            netDict[i[1]] = set()
        
        netDict[i[0]].add(i[1])
        netDict[i[1]].add(i[0])
    
    print("writing")
    count = 0
    length = len(geneList)
   
    #mat = np.zeros((len(geneList),len(geneList)))
    mat = []
   
    for i in range(len(geneList)):
        if count% 5000 == 0:
            print ("...%i/%i"%(count,length))
        count +=1
        mat.append([])
        for k in range(len(geneList)):
            if geneList[k] in netDict[geneList[i]]:
                mat[-1].append(1) 
                #print("found")
            else:
                mat[-1].append(0) 
    mat = np.array(mat)
    return mat,geneList
    
def readGWAS(gwasPath geneList, geneCol =0, pvalCol = 1, sep = "\t", header = True)):

    pvalDict = {}
    with open(gwasPath) as data:
        for i in data.readlines():
            if header:
                header = False
                continue
            i = i.strip().split(sep)
            pvalDict[i[geneCol]] = float(i[pvalCol])
    
    heats =[]
    heatTups = []
    for i in geneList:
        if i not in pvalDict:
            heats.append(0)
        else:
            heats.append(-math.log(pvalDict[i],1.5))
        heatTups.append((heats[i],geneList[i])
    
    heatTups.sort(reverse =True)
    top25 = [x[1] for x in heatTups[:25]
    
    return heats,top25


def heatProp (adjMatrix, heat, prop = 3, time = 1):
    degree = np.sum(adjMatrix,1)
    heatProp = np.identity(adjMatrix.shape[0])
    laplace = np.identity(adjMatrix.shape[0])
    
    print("Calculating Laplace")
    for i in range(len(degree)):
        laplace[i,i] = degree[i]
    
    laplace = adjMatrix - laplace
    laplace = laplace*time
    
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
    
def output_Arya(heatRanks, output):
    with open(output, 'w') as fileout:
        for tups in heatRanks:
            gene = tup[0]
            heat = tup[1]
            fileout.write("%s\t%s" % (gene, heat))
            
def output(heat,geneList,top25, output):
    with open(output, 'w') as fileout:
        for g in top25:
            fileout.write("%s\n"%g)
            
        heatTups = []
        for i in range(len(heat)):
            heatTups.append((heat[i],geneList[i]))
            
        heatTups.sort(reverse=True)
        count = 0
        for h in heatTups:
            if h[1] not in top25:
                fileout.write("%s\n"%h[1])
                count +=1
            if count >= 975:
                break
             
        
def maxheats(masterheatRanks):

    for h in masterheatRanks:
            #curr = [x.strip() for x in open(f).readlines()]
            if len(maxList) ==0:
                for i in h:
                    maxList.append([i])
            else:
                for i in range(len(h)):
                    maxList[i].append(curr[i])
                    
    finalheats = []
    for x in maxList:
        finalheats.append(max(x))
    return finalheats
def getTimes (maxT):

    times = [.01,.1,.25,.5,.75,1,3,5,10,.001,.0001]
    priorityTups = [(5,0),(1,1),(8,2),(6,3),(9,4),(2,5),(7,6),(3,7),(4,8),(10,9),(11,10)]
    priorityTups.sort()
    
    finalTimes = []
    for i in priorityTups:
        if len(finalTimes) >= 5:
            break
        if times[i[1]] <= maxT:
            finalTimes.append(times[i[1]])
    return finalTimes
    
def main():
    args = parseargs(sys.argv[1:])
    adjMatrix,geneList = readNetwork(args.network)
    heat,top25 = readGWAS(args.gwas,geneList)
    timeVals = getTimes(args.uppertimebound)
    masterheatranks = []
    for t in timeVals:
        heatranks = heatProp(adjMatrix, heat,time = t)
        masterheatranks.append(heatranks)
    finalheats = maxheats(masterheatranks)
    output(finalheats,geneList,top25, args.output)

if __name__=="__main__":
	main()
