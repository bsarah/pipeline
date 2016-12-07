
import subprocess
from os import listdir

def endSlash(tupleOfArgs):
    '''
    takes a tuple of Args. 
    returns a tuple of args with '/' at the end
    '''

    returnlist = list()
    for arg in tupleOfArgs:
        if not arg.endswith('/'):
            returnlist.append(arg+'/')
        else:
            returnlist.append(arg)
    return tuple(returnlist)

def handlePathToRepo(pathToRepo):
    if not pathToRepo.endswith('/pipeline/'):
        if pathToRepo.endswith('/scripts_cam/'):
            return pathToRepo
        else:
            return pathToRepo + 'pipeline/scripts_cam/'
    else:
        return pathToRepo + 'scripts_cam/'

def makeSpeciesList(_list, genomes, reference):
    for f in listdir(genomes):
        if isfile(join(genomes,f)):
            species = f.split('.')[0]
            if species != reference:
                _list.append(species)
    return _list
                
def doesOutputDirExist(outPutDir):
    try:
        listdir(outPutDir)
    except FileNotFoundError:
        subprocess.call('mkdir '+outPutDir, shell=True)

if __name__ == '__main__':
    pass
