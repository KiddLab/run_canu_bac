import sys
import os
import gzip
import shutil
import subprocess
import datetime
import numpy as np




###############################################################################
def add_breaks_to_line(seq,n=50):
    myList = []
    myList = [i for i in seq]
    newList = []
    c = 0
    for i in myList:
        newList.append(i)
        c += 1
        if c % n == 0 and c != (len(myList)):
            newList.append('\n')
    myStr = ''.join(newList)
    return myStr    
###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed', flush=True)
        print(cmd)
        sys.exit(1)
###############################################################################
# Helper function to run commands,
# doesn't check return or print to log.  Use for 'grep' so that doesn't
# fail if no items found
def runCMDNoFail(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
###############################################################################
def get_4l_record(myFile):
    #fastq style file...
    # just return sequence len
    myLine1 = myFile.readline()
    if myLine1 == '':
        return ''
    myLine2 = myFile.readline()
    myLine3 = myFile.readline()
    myLine4 = myFile.readline()
    return [myLine1,myLine2,myLine3,myLine4]
###############################################################################
# setup paths to default programs to use and checks for required programs
def check_prog_paths(myData):    
    print('Checking samtools...')
    if shutil.which('samtools') is None:
        print('samtools not found in path! please fix (module load?)', flush=True)
        sys.exit()

    print('Checking blat...')
    if shutil.which('blat') is None:
        print('blat not found in path! please fix (module load?)', flush=True)
        sys.exit()

    print('Checking minimap2...')
    if shutil.which('minimap2') is None:
        print('minimap2 not found in path! please fix (module load?)', flush=True)
        sys.exit()

    print('Checking canu...')
    if shutil.which('canu') is None:
        print('minimap2 not found in path! please fix (module load?)', flush=True)
        sys.exit()



#####################################################################
def parse_paf_line(line):
    pafLine = {}
    pafLine['qName'] = line[0]
    pafLine['qLen'] = int(line[1])
    pafLine['qStart'] = int(line[2])
    pafLine['qEnd'] = int(line[3])
    pafLine['strand'] = line[4]
    pafLine['tName'] = line[5]
    pafLine['tLen'] = int(line[6])
    pafLine['tStart'] = int(line[7])
    pafLine['tEnd'] = int(line[8])
    pafLine['numMatch'] = int(line[9])
    pafLine['alignBlockLen'] = int(line[10])
    pafLine['mapQ'] = int(line[11])
    return pafLine

##############################################################################
# write initial info to log file
def write_initial_log(myData):
    now = datetime.datetime.now()
    myData['logFile'].write(str(now))
    myData['logFile'].write('\n')
    myData['logFile'].write('Initial config information\n')    
    k = myData.keys()
    k = list(k)
    k.sort()
    toSkip = ['logFile']
    for i in k:
        if i in toSkip:
            continue
        myData['logFile'].write('%s\t%s\n' % (i,myData[i]))
    myData['logFile'].write('\n')        
#####################################################################
def filter_contam_longread(myData):
    myData['longReadContam'] = myData['outDirBase'] + 'longread.contam.paf'
    myData['longReadFilt'] = myData['outDirBase'] + 'longread.filt.fq.gz'   
    myData['longReadFiltFail'] = myData['outDirBase'] + 'longread.fail.fq.gz'   
    
    # check to see if already ran
    if os.path.isfile(myData['longReadFilt']) and os.path.getsize(myData['longReadFilt']) > 0:
        print('Already ran long read filter contamination!', flush=True)
        print('Will continue anyway')
    
    if myData['longreadtype'] == 'ont':
       mapX = 'map-ont'
    else:
       print('uknown long read type, option not clear', flush=True)
       print(myData['longreadtype'])
       sys.exit()
    
    cmd = 'minimap2 -t 1 -x %s %s %s > %s ' % (mapX,myData['contam'],myData['longread'],myData['longReadContam'] )
    print(cmd)
    runCMD(cmd)
    
    toDrop = {}
    inFile = open(myData['longReadContam'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)
        
        if((pafLine['numMatch']/pafLine['qLen']) > 0.25):
            toDrop[pafLine['qName']] = 0
    inFile.close()
    print('Read in %i names of long reads' % len(toDrop), flush=True)
    
    inFile = gzip.open(myData['longread'],'rt')
    outPass = gzip.open(myData['longReadFilt'],'wt')
    outFail = gzip.open(myData['longReadFiltFail'],'wt')
    
    numDrop = 0
    numWrite = 0
    
    while True:
        r1 = get_4l_record(inFile)
        if r1 == '':
            break
        
        n1 = r1[0].rstrip()
        n1 = n1[1:].split()[0]
        
        if n1 in toDrop:
            numDrop +=1
            outFail.write('%s%s%s%s' % (r1[0],r1[1],r1[2],r1[3]))

        else:
            numWrite +=1
            outPass.write('%s%s%s%s' % (r1[0],r1[1],r1[2],r1[3]))
    inFile.close()
    outPass.close()
    outFail.close()
    
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nSearched for contamination in long reads!\n')
        fstream.write('Total long reads:\t%i\n' % (numDrop+numWrite) )
        fstream.write('Removed\t%i\t%f\n' % (numDrop, numDrop/(numDrop+numWrite)))
        fstream.write('Kept\t%i\t%f\n' % (numWrite, numWrite/(numDrop+numWrite)))
        fstream.flush()
#####################################################################
# returns N50, N90, Nk etc, assumes that data is reverse sorted, k is integer
# so will do k/100
def get_nx(data,totBp,k):
    targetBp = (k/100.0) * totBp
    sum = 0
    for i in range(0,len(data)):
        sum += data[i]
        if sum >= targetBp:
            break
    return data[i]
#####################################################################
def print_longread_stats(fqFileName,myData):
    readLens = []
    inFile = gzip.open(fqFileName,'rt')
    while True:
        rec = get_4l_record(inFile)
        if rec == '':
            break
        seq = rec[1].rstrip()
        readLens.append(len(seq))
                    
    inFile.close()
    numSeqs = len(readLens)
    totLen = sum(readLens)
    readLens.sort(reverse=True)
    med = np.median(readLens)
    avg = np.mean(readLens)
    
    
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nstats for reads in\t%s\n' % fqFileName)
        fstream.write('total reads\t%i\n' % (numSeqs))
        fstream.write('total bp\t%i\n' % (totLen))
        fstream.write('median\t%i\n' % (med))
        fstream.write('mean\t%f\n' % (avg))
        fstream.write('n50\t%i\n' % get_nx(readLens,totLen,50))
        fstream.write('5 longest read lengths:\n')
        for i in range(5):
            fstream.write('\t%i\n' % readLens[i])        
        fstream.flush()        
#######################################################################
def run_canu_assem(myData):  
    myData['canuDir'] = myData['outDirBase'] + 'canuDir'
    cmd = 'canu -useGrid=false -maxThreads=%i -corThreads=1 ' % (myData['numThreads'])
    cmd += ' -d %s ' % myData['canuDir']
    cmd += ' -p %s ' % myData['name']
    cmd += ' genomeSize=250k '
    cmd += ' -nanopore %s ' % myData['longReadFilt']
    
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nCanu cmd is:\n%s\n' % cmd)
        fstream.flush()
        
    runCMD(cmd)    
#######################################################################    









