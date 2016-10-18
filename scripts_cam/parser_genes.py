import subprocess
from os import listdir
from os.path import isfile, join
from dataStructures_lessMem import Gene

'''
files given to the gene parser should be: 
  -in the modified BED format below.
  -end with .bed (or .bed.gz, .bed.Z, .bed.bz2)
  -zipped in gzip or unzipped
  -species in list must match those in the MAF files given

File Format:
  -Heading (optional):
    -line must start with a '#'
  -Body:
    -one line = one gene
    -tab separated elements
    -chromosome     species     startCoord     endCoord     +or-     sequence     secondary_structure(optional)
'''



def parseGenes(filePath, outputDir, speciesList=list()):
    listOfGenes = list()
    geneNumber = 0

    for fileName in listdir(filePath):
        filePlusPath= join(filePath, fileName)
        if isfile(filePlusPath):
            f = catchExceptions(filePlusPath, outputDir)

            for line in f:
                if line[0] in '# \n':
                    pass
                else:
                    geneNumber += 1
                    lineList = line.split()
                    if lineList[1] not in speciesList:
                        speciesList.append(lineList[1])
                    if len(lineList) == 7:
                        structure = lineList[6]
                    else:
                        structure = None
                        #                       panTro       chr3         84545             85
                        listOfGenes.append(Gene(lineList[1], lineList[0], int(lineList[2]), getGeneLength(int(lineList[2]), int(lineList[3])),\
                                                lineList[4], geneNumber, structure, lineList[5]))
                        #                       +            32          ((___>><<) ATTCGTAGCAT

    return listOfGenes, speciesList



def catchExceptions(fileName, outputDir):
    '''
    checks if file is in correct format( '.bed', '.bed.gz', '.bed.Z', '.bed.bz2')
    checks if file is zipped
    checks if file can be opened
    checks if file can be read

    unzips file if zipped
    '''
    
    name = fileName.split('/')[-1]#fileName without pathway ex: test1.bed.gz
    
    #proper extension
    if fileName.endswith('.bed'):
        pass

    elif fileName.endswith(('.bed.gz','.bed.Z','.bed.bz2')):

        unzippedDir = outputDir+'unzippedGenes/'
        if 'unzippedGenes' in listdir(outputDir):
            subprocess.call("rm -r "+unzippedDir, shell=True)
        subprocess.call('mkdir '+unzippedDir,shell=True)
        
        #unzippedFileName = ('/'.join(fileName.split('/')[:-1]))+'/'+name.split('.')[0]+'.maf'
        unzippedFileName = unzippedDir + name.split('.')[0]+'.bed'
        print("unzipping maf file...")
        subprocess.call('/usr/bin/zcat '+fileName+' > '+unzippedFileName, shell=True)
        print("done")
        fileName = unzippedFileName
    else:
        raise Exception("{} not of bed format. Use a .bed file or gziped .bed file('.bed.gz','.bed.Z','.bed.bz2')".format(fileName))

    #can unzipped file be opened
    #try:
    f = open(fileName, 'r')
    #except IOError:
        #raise Exception("cannot open {}".format(fileName))
    #else:
    return f

def getGeneLength(start, end):
    '''
    returns the length of the gene
    
    if we are on the + strand we read left to right so end is bigger so subtract start from end.
    If we are on the - strand we read right to left so start is bigger so subtract end from start.
    '''
    if end > start:
        return end - start
    return start - end


if __name__ == '__main__':

    geneObjects = parseGenes('../textFiles/test/','../textFiles/outPut')
    print("len: {}".format(len(geneObjects)))
    print(geneObjects)