'''
This is the shell script that will run the gene evolution pipeline.

It will take the following arguments:

  1)File pathway to the Genomes that are to be compared
  2)File pathway to the Multiple sequence allignments of those genomes
      These allignments must be labeled identically to the genomes (ex: genomename: panTro, multiple sequence alignment name: panTro)
      The allignment format should be:
        -MAF (file extension '.maf')
          -uncompressed OR
          -compressed with gzip (file extension '.maf.gz', '.maf.Z', '.maf.bz2')

  3)File pathway to a directory where you would like the output to be placed
  4)Name of the reference species used in the multiple sequence alignment as it appears in the multiple sequence alignment / genomes (these should be the same)

Optional Arguments:



  -sg:Genes used in analysis will be found by cmsearch. User must provide a model of the genes they are looking for
  -og:Genes used in analysis are those provided by the user. There should be one file for every species. The files should be named:
          'panTro.bed' or 'dm6.bed', etc.
        Files should be in bed format, with one gene per line of the file.
  -q :Remove a percentage of the genes 

It will perform the following steps:

  1)
'''
import argparse #python module dealing with command line arguments
import sys
import subprocess #python module allowing other non-python files to be called
from os import listdir #listdir returns a list of everything in a directory
from os.path import isfile, join #join joins a file extension with a file name
from parser_MAF_lessMem import maf2bed #maf2bed parses maf files and writes bed files
from parser_infernal import parseInfernal #parses infernal files and returns a list of Genes
from parser_genes import parseGenes #parses user provied gene files and returns a list of Genes


#globals
#locationInfernal = '/opt/bin'
#locationModel = '/scr/k61san/trnaevo/CMs/trna.cm'

parser = argparse.ArgumentParser()
parser.add_argument("genomes", type=str, help="directory where the genomes are")
parser.add_argument("multiSeq", type=str, help="directory where the multiple sequence allignments are")
parser.add_argument("infernalPath", type=str, help="path to infernal")
parser.add_argument("outPutDir", type=str, help="directory where you want the output files to be placed")
parser.add_argument("referenceSpecies", type=str, help="the name of the reference species as it appears in the "\
                    + "multiple sequence alignment.")
parser.add_argument("-sg","--search_genes", action="store", help="search for genes using infernal cmsearch. "\
                    + "File pathway provided must lead to a rna model.")
parser.add_argument("-og","--own_genes", action="store", help="use prespecified genes in analysis. File path "\
                    + "provided should be a directory. Files in this directory should be named for the species "\
                    + "that it concerns. The names of the species should match those in the multiple sequence "\
                    + "alignments exactly. The files should contain the genes of interest in BED format.")
parser.add_argument("-q","--quality", action="store", default=0, type=int, help="a number "\
                    + "specifying the percentage of blocks that will be thrown away. These blocks will be the lowest "\
                    + "quality blocks. Warning: throwing away too many blocks may lead to an inconclusive analysis.")
parser.add_argument("-incE", action="store", type=float, help="consider sequences <= this E-value threshold as "\
                    + "significant")
parser.add_argument("-incT", action="store", type=float, help="consider sequences >= this score threshold as significant")

args = parser.parse_args()

if args.infernalPath.endswith('/'):
    locationInfernal = args.infernalPath.rstrip('/')
else:
    locationInfernal = args.infernalPath
    
args.quality = int(args.quality)
if args.quality < 0 or args.quality > 99:
    raise Exception("quality value must be a number between 0 and 99 corresponding to the percentage of "\
                    "blocks that will be thrown away.")

if args.own_genes != None and args.search_genes != None:
    raise Exception("both own_genes and search_genes given. own_genes and search_genes are mutually exclusize.")

if not args.outPutDir.endswith('/'):
    args.outPutDir = args.outPutDir + '/'


#returns a list of all the files in a given directory with the full file path
genomeFiles = [join(args.genomes,f) for f in listdir(args.genomes) if isfile(join(args.genomes,f))]
multiSeqFiles = [join(args.multiSeq,f) for f in listdir(args.multiSeq) if isfile(join(args.multiSeq,f))]

listOfSpecies = [args.referenceSpecies]
for f in listdir(args.genomes):
    if isfile(join(args.genomes,f)):
        species = f.split('.')[0]
        if species != args.referenceSpecies:
            listOfSpecies.append(species)

#listOfSpecies f.split('.')[0] for f in listdir(args.genomes) if isfile(join(args.genomes,f))]

# ---------------------------------THE FIRST SPECIES IN THE LIST MUST BE THE REFERENCE SPECIES----------------------------------------------
#listOfSpecies = ['dm6', 'anoGam1', 'apiMel4', 'droAlb1', 'droAna3', 'droBia2', 'droBip2', 'droEle2', 'droEre2', 'droEug2',\
#                 'droFic2', 'droGri2', 'droKik2', 'droMir2', 'droMoj3', 'droPer1', 'droPse3', 'droRho2', 'droSec1', 'droSim1',\
#                 'droSuz1', 'droTak2', 'droVir3', 'droWil2', 'droYak3', 'musDom2', 'triCas2'] #------------------------Testing -----------------------------

if len(multiSeqFiles) == 0:
    raise Exception("No multiple sequence alignments given")

try:
    listdir(args.outPutDir)

except FileNotFoundError:
    subprocess.call('mkdir '+args.outPutDir, shell=True)
    
if args.search_genes != None:#if we need to search for genes
    
    #prep
    if len(genomeFiles) == 0:
        raise Exception("No genomes given. Make sure all paths given are directories not a files")
    if 'infernalIn' in listdir(args.outPutDir):
        subprocess.call("rm -r "+args.outPutDir+'infernalIn', shell=True)
    subprocess.call("mkdir "+args.outPutDir+'infernalIn', shell=True)

    i = 1 #counter for the number of genomes searched

    if args.incE == None:
        if args.incT == None:
            infernalOptArgs = ''
        else:
            infernalOptArgs = ' --incT '+str(args.incT)
    else:
        if args.incT == None:
            infernalOptArgs = ' --incE '+str(args.incE)
        else:
            infernalOptArgs = ' --incE '+str(args.incE)+' --incT '+str(args.incT)
    #print("incE: {}".format(str(args.incE)))
    #print("incT: {}".format(str(args.incT)))
    #print("optArgs: {}".format(infernalOptArgs))
            
    #search the genomes for genes
    for g_file in genomeFiles:

        species = g_file.split('/')[-1].split('.')[0]
        print("searching {} for genes of interest...".format(g_file))
        subprocess.call(locationInfernal+'/cmsearch -o '+args.outPutDir+'infernalIn/'+species+'.out '+infernalOptArgs+' '+args.search_genes+' '+g_file, shell=True)
                
        print("done searching genome {}/{}".format(i,len(genomeFiles)))
        i+=1

    print("converting infernal files to list of gene objects...")
    geneObjects, versionInfo, listOfSpecies = parseInfernal(args.outPutDir+'/infernalIn', [args.referenceSpecies])
    print(listOfSpecies)
    print("done.")

elif args.own_genes != None:
    print("parsing given genes...")
    geneObjects, listOfSpecies = parseGenes(args.own_genes,args.outPutDir, [args.referenceSpecies])
    versionInfo = "User provided genes"
    print("done")

else:
    raise Exception("either own_genes or search_genes must be given")


if len(multiSeqFiles) == 0:
    raise Exception("No multiple sequence allignments given")

print("sorting multiple sequence alignment files...")
multiSeqFiles.sort()#makes the output files alphabetical
print("done")
print("parsing maf files, writing valid blocks in bed format...")

maf2bed(multiSeqFiles, args.outPutDir, geneObjects, listOfSpecies, args.quality, versionInfo)
#print("done. Bed files stored in {}".format(('/'.join(multiSeqFiles[0].split('/')[:-2]))+'/bed'))

