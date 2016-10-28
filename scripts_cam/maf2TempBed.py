import sys
import subprocess
import argparse
from os import listdir
from os.path import isfile, join
from dataStructure_lessMem import catchExceptions

parser = argparse.ArgumentParser()
parser.add_argument("blockNum", type=int, help="current block number between all maf files")
parser.add_argument("mafDir", type=str)
parser.add_argument("iteration", type=int)
parser.add_argument("writeFilesDir", type=str)

args = parser.parse_args()

blockNum = args.blockNum
iteration = args.iteration
listSpeciesFiles = [join(args.writeFilesDir, f) for f in listdir(args.writeFilesDir) if isfile(join(args.writeFilesDir, f))]
listMafFiles = [join(args.mafDir, f) for f in listdir(args.mafDir) if isfile(join(args.mafDir, f))]

for line in sys.stdin:
    
    print("converting maf files to temp files...")
    lastLineA = False

    if line[0] not in {'#','a','s','i','e','q',' ','\n'}:
        print("Ignoring {}. File not of maf format".format(readFile.name))
        break

    if line[0] == 'a':
        blockNum += 1
        lastLineA = True
        buf = line.split()[1]
               #read score, with and without a decimal place
                                                                       
        if '.' in buf.split('=')[1]:
            score = int(float(buf.split('=')[1]))
        else:
            score = int(buf.split('=')[1])

    if line[0] == 's':
        buf = line.split()
        ID = buf[1].split('.')

        try:
            test = ID[1]
        except IndexError:
            print("Ignoring {}. File lacks chromosome information".format(readFile.name))
            break #break out of for line in readFile for loop, continue to next filename in mafdir (recursive call)                                                                         
        else:
            i = 0
            for writeFile in listSpeciesFiles:
                writeFile = open(writeFile, 'w')
                i+=1
                numSearches = len(listSpeciesFiles)

                writeFileName = writeFile.name.split('/')[-1]
                if writeFileName.startswith(ID[0]):#ex: does dm6_temp.bed start with dm6                                                                                                            
                    writeFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                                    ID[1],ID[0]+'_'+str(blockNum), buf[2], buf[3], buf[4], lastLineA, score))
                    #               chr   species_blocknum         start   length  strand  isRefSpecies allignmentScore

                elif i > numSearches:
                    print("writeFileName: {}".format(writeFileName))
                    raise Exception("the species of a maf allignment was not in the list of species given\n"+\
                                    "allignment species: {}".format(ID[0]))
                lastLineA = False
    print("maf file read")

#Call next iteration
#1 Recursive call for each maf file
catchExceptions(listMafFiles[args.iteration+1])
subprocess.call("zcat "+listMafFiles[args.iteration+1]+" | python3 maf2TempBed.py "+str(blockNum)+" "+args.mafDir+" "+str(iteration+1))
