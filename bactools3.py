import sys
import os
import gzip
import shutil
import subprocess
import datetime
import numpy as np
import matplotlib.pyplot as plt




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
    pafLine['tags'] = line[12:]    
    return pafLine

##############################################################################
def paf_line_is_primary(pafLine):
    for tag in pafLine['tags']:
        if tag[0:3] == 'tp:':
            if tag.split(':')[2] == 'P':
                return True
            else:
                return False
    print('No tp tag found, unclear!')
    sys.exit()
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
        return
    
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
def run_miropeats_self(myData):
   # run miropeats intra 

    print('Checking miropeats...')
    if shutil.which('miropeats') is None:
        print('miropeats not found in path! please fix (module load?)', flush=True)
        sys.exit()

    if shutil.which('ps2pdf') is None:
        print('ps2pdf not found in path! please fix (module load?)', flush=True)
        sys.exit()

    myData['miroOut'] = myData['outDir'] + myData['name'] + '.mirropeats.200.sel.out'
    myData['miroOutPS'] = myData['outDir'] + myData['name'] + '.mirropeats.200.sel.out.ps'
    myData['miroOutPDF'] = myData['outDir'] + myData['name'] + '.mirropeats.200.sel.out.pdf'    
    
    
    cmd = 'miropeats -s 200 -onlyintra -o %s -seq %s > %s ' % (myData['miroOutPS'],myData['contig'],myData['miroOut'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nmiropeats cmd is:\n%s\n' % cmd)
        fstream.flush()
        
    runCMD(cmd)        
    cmd = 'ps2pdf %s %s' % ( myData['miroOutPS'],myData['miroOutPDF'])
    runCMD(cmd)        
#######################################################################    
def map_to_contig_paf(myData):
    if shutil.which('minimap2') is None:
        print('miropeats not found in path! please fix (module load?)', flush=True)
        sys.exit()
        
    myData['PAFout'] = myData['outDir'] + 'align.paf'
    myData['PAFbed'] = myData['outDir'] + 'align.paf.bed'    
    
    
    cmd = 'minimap2 -x '
    if myData['longreadtype'] == 'ont':
        cmd += ' map-ont'
    else:
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nERRROR!!!\nlong read type is :\n%s\nI do not know it' % myData['longreadtype'])            
            fstream.flush()
        sys.exit() # fail!
    
    cmd += ' -c %s %s > %s ' % (myData['contig'],myData['longread'],myData['PAFout'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\minimap2 cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)        
    
    outFile = open(myData['PAFbed'],'w')
    inFile = open(myData['PAFout'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)        
        if paf_line_is_primary(pafLine) is True and pafLine['mapQ'] != 0:
            nl = [pafLine['tName'],pafLine['tStart'],pafLine['tEnd'],pafLine['qName']]
            nl = [str(j) for j in nl]
            nl = '\t'.join(nl) + '\n'
            outFile.write(nl)
    outFile.close()
    inFile.close()
#######################################################################    
def make_windows_bed(myData,ws,stepsize):
    
    # first, need the .fai file
    faSeqs = read_fasta_file_to_dist(myData['contig'])
    
    if len(faSeqs) != 1:
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nERROR\nHave to many seqs in %s\n' % myData['contig'])
            fstream.flush()
        sys.exit()
    
    
    myData['contigLenFile'] = myData['outDir'] + 'contig_len.txt'
    outFile = open(myData['contigLenFile'],'w')
    for seq in faSeqs:
        outFile.write('%s\t%i\n' % (seq,faSeqs[seq]['seqLen']))
        myData['contigLen'] = faSeqs[seq]['seqLen']
    outFile.close()
    
    # now make the windows
    myData['bedWindowsFile'] = myData['outDir'] + myData['name'] + '.%i.%i.bed' % (ws,stepsize)
    
    cmd = 'bedtools makewindows -g %s -w %i -s %i > %s' % (myData['contigLenFile'], ws,stepsize,myData['bedWindowsFile'])
    runCMD(cmd)
#######################################################################    
def read_fasta_file_to_dist(fastaFile):
    myDict = {}
    inFile = open(fastaFile,'r')
    line = inFile.readline()
    line = line.rstrip()
    if line[0] != '>':
        print('ERROR, FILE DOESNNOT START WITH >', flush=True)
        sys.exit()
    myName = line[1:].split()[0]
    myDict[myName] = {}
    myDict[myName]['seq'] = ''
    myDict[myName]['seqLen'] = 0    
    mySeq = ''
    while True:
        line = inFile.readline()
        if line == '':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            break
        line = line.rstrip()
        if line[0] == '>':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            myName = line[1:].split()[0]
            myDict[myName] = {}
            myDict[myName]['seq'] = ''
            myDict[myName]['seqLen'] = 0    
            mySeq = ''
            continue
        mySeq += line
    inFile.close()
    return myDict
###############################################################################
def make_coverage_plot(myData):
    # get the mean coverage in windows
    myData['coverageInWindows'] = myData['bedWindowsFile'] + '.coverage'
    myData['coverageInWindowsPlt'] = myData['bedWindowsFile'] + '.coverage.png'
    
    
    cmd  = 'bedtools coverage -mean  -a %s -b %s > %s ' % (myData['bedWindowsFile'],myData['PAFbed'],myData['coverageInWindows'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\ncoverage cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)

    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\ncontig len is:\t%i\n' % myData['contigLen'])
        fstream.flush()

    
    xpos = []
    depths = []
    inFile = open(myData['coverageInWindows'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        b = int(line[1])
        e = int(line[2])
        meanDepth = float(line[3])
        mp = (e-b)/2 + b
        xpos.append(mp)
        depths.append(meanDepth)
   
    inFile.close()
        
    plt.figure(figsize=(8,6))        
    plt.plot(xpos,depths,'o',color='black',markersize=3)
    plt.xlim([1,myData['contigLen']])    
    plt.xlabel(myData['name'])
    plt.ylabel('Read Depth')    
    plt.savefig(myData['coverageInWindowsPlt'])
    plt.close()
#######################################################################    
def make_gc_plot(myData):
    # get the mean coverage in windows
    myData['GCInWindows'] = myData['bedWindowsFile'] + '.gc'
    myData['GCInWindowsPlt'] = myData['GCInWindows'] + '.png'
    
    cmd = 'bedtools nuc -fi %s -bed %s > %s ' % (myData['contig'],myData['bedWindowsFile'],myData['GCInWindows'])
    
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\ngc cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)
    
    xpos = []
    gc = []
    inFile = open(myData['GCInWindows'],'r')
    for line in inFile:
        if line[0] == '#':
            continue            
        line = line.rstrip()
        line = line.split()
        b = int(line[1])
        e = int(line[2])
        gcW = float(line[3])
        mp = (e-b)/2 + b
        xpos.append(mp)
        gc.append(gcW)   
    inFile.close()
        
    plt.figure(figsize=(8,6))        
    plt.plot(xpos,gc,'o',color='blue',markersize=3)
    plt.xlim([1,myData['contigLen']])  
    plt.ylim([0,1.0])  
    plt.xlabel(myData['name'])
    plt.ylabel('GC Fraction')    
    plt.savefig(myData['GCInWindowsPlt'])
    plt.close()    
#######################################################################    
def make_coverage_plot_showlong(myData,numToShow):
    # get the mean coverage in windows
    myData['coverageInWindowsShowLongPlt'] = myData['bedWindowsFile'] + '.coverage.showlong.png'
    
    
    
    xpos = []
    depths = []
    inFile = open(myData['coverageInWindows'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        b = int(line[1])
        e = int(line[2])
        meanDepth = float(line[3])
        mp = (e-b)/2 + b
        xpos.append(mp)
        depths.append(meanDepth)   
    inFile.close()
    
    # need to read in the segments of reads.
    alignSegs = []
    inFile = open(myData['PAFbed'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        b = int(line[1])
        e = int(line[2])
        alignSegs.append([b,e])    
    inFile.close()
    
    # sort be length
    alignSegs.sort(key=lambda i: i[1]-i[0], reverse=True)
    print(alignSegs[0:numToShow])
    
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\naligned position of %i longest reads\n' % (numToShow) )
        for i in range(numToShow):
            fstream.write('%i\t%i\t%i\n' % (alignSegs[i][0],alignSegs[i][1],alignSegs[i][1]-alignSegs[i][0] ) )
        fstream.flush()
        
    plt.figure(figsize=(8,6))        
    plt.plot(xpos,depths,'o',color='black',markersize=3)
    plt.xlim([1,myData['contigLen']])    
    plt.xlabel(myData['name'])
    plt.ylabel('Read Depth')   
    plt.title('Coverage with %i Longest Aligned Segments' % numToShow) 
    
    for i in range(numToShow):
        b = alignSegs[i][0]
        e = alignSegs[i][1]
        y = i * 1.25 + 1
        plt.plot([b+1,e],[y,y],color='red')
    
    
    plt.savefig(myData['coverageInWindowsShowLongPlt'])
    plt.close()
#######################################################################    










