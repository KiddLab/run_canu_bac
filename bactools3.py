import sys
import os
import gzip
import shutil
import subprocess
import datetime
import numpy as np
import matplotlib.pyplot as plt


##############################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c   
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''    
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
##############################################################################
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
        print('canu not found in path! please fix (module load?)', flush=True)
        sys.exit()

    print('Checking miropeats...')
    if shutil.which('miropeats') is None:
        print('miropeats not found in path! please fix (module load?)', flush=True)
        sys.exit()

    print('Checking ps2pdf...')
    if shutil.which('ps2pdf') is None:
        print('ps2pdf not found in path! please fix (module load?)', flush=True)
        sys.exit()

    print('Checking racon...')
    if shutil.which('racon') is None:
        print('racon not found in path! please fix (module load?)', flush=True)
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
    myData['miroOut'] = myData['outDir'] + myData['name'] + '.mirropeats.200.sel.out'
    myData['miroOutPS'] = myData['outDir'] + myData['name'] + '.mirropeats.200.sel.out.ps'
    myData['miroOutPDF'] = myData['outDir'] + myData['name'] + '.mirropeats.200.sel.out.pdf'    

    if os.path.isfile(myData['miroOutPS'])  is True:
        cmd = 'rm %s'  % myData['miroOutPS']
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nprevious miroOutPS file found, removing it\n%s\n' % cmd)
            fstream.flush()
        runCMD(cmd)

    if os.path.isfile(myData['miroOutPDF'])  is True:
        cmd = 'rm %s'  % myData['miroOutPDF']
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nprevious miroOutPDF file found, removing it\n%s\n' % cmd)
            fstream.flush()
        runCMD(cmd)

    
    cmd = 'miropeats -s 200 -onlyintra -o %s -seq %s > %s ' % (myData['miroOutPS'],myData['contig'],myData['miroOut'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nmiropeats cmd is:\n%s\n' % cmd)
        fstream.flush()
        
    runCMD(cmd)        
    cmd = 'ps2pdf %s %s' % ( myData['miroOutPS'],myData['miroOutPDF'])
    runCMD(cmd)        
#######################################################################    
def map_to_contig_paf(myData):        
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
        fstream.write('\nminimap2 cmd is:\n%s\n' % cmd)
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
    faSeqs = read_fasta_file_to_dict(myData['contig'])
    
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
def read_fasta_file_to_dict(fastaFile):
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
    
    maxDepth = 1.1 * max(depths) # 10% more
    
        
    plt.figure(figsize=(8,6))        
    plt.plot(xpos,depths,'o',color='black',markersize=3)
    plt.xlim([1,myData['contigLen']])    
    plt.ylim([0,maxDepth])
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
    
    
    maxDepth = 1.1 * max(depths) # 10% more

        
    plt.figure(figsize=(8,6))        
    plt.plot(xpos,depths,'o',color='black',markersize=3)
    plt.xlim([1,myData['contigLen']])  
    plt.ylim([0,maxDepth])
      
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
def make_vector_esp_plot(myData):
    print('ready to do vector and esp')
    
    myData['vectorMapOut'] = myData['outDir'] + 'vector-map.paf'
    
    if os.path.isfile(myData['vectorMapOut'])  is True:
        cmd = 'rm %s'  % myData['vectorMapOut']
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nprevious vector map out file found, removing it\n%s\n' % cmd)
            fstream.flush()
        runCMD(cmd)
    
    cmd = 'minimap2 -c -x asm5 %s %s > %s ' % (myData['contig'],myData['vector'],myData['vectorMapOut'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nmap vector cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)
    
    # read in vector hits
    vectorHits = []
    inFile = open(myData['vectorMapOut'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)   
        vectorHits.append([pafLine['tStart'],pafLine['tEnd'],pafLine['strand']])        
    inFile.close()
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('Vector hits:\n')
        for v in vectorHits:
            fstream.write('%s\t%s\t%s\n' % (v[0],v[1],v[2]))
        fstream.flush()
    
    # get R1 and R2 hits
    myData['R1fa'] = myData['outDir'] + 'R1.fa'
    myData['R1Map'] = myData['outDir'] + 'R1-map.paf'

    myData['R2fa'] = myData['outDir'] + 'R2.fa'
    myData['R2Map'] = myData['outDir'] + 'R2-map.paf'

    if os.path.isfile(myData['R1Map'])  is True:
        cmd = 'rm %s' % myData['R1Map']
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nprevious R1Map map out file found, removing it\n%s\n' % cmd)
            fstream.flush()
        runCMD(cmd)

    if os.path.isfile(myData['R2Map'])  is True:
        cmd = 'rm %s' % myData['R2Map']
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nprevious R2Map map out file found, removing it\n%s\n' % cmd)
            fstream.flush()
        runCMD(cmd)


    
    cmd = 'samtools faidx %s %s > %s' % (myData['libraryEndSeqsFA'],myData['cloneName'] + '_R1', myData['R1fa'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nextract and map R1 and R2:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)
    cmd = 'samtools faidx %s %s > %s' % (myData['libraryEndSeqsFA'],myData['cloneName'] + '_R2', myData['R2fa'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)

    cmd = 'minimap2 -c -x asm5 %s %s > %s ' % (myData['contig'],myData['R1fa'],myData['R1Map'])
    runCMD(cmd)
    cmd = 'minimap2 -c -x asm5 %s %s > %s ' % (myData['contig'],myData['R2fa'],myData['R2Map'])
    runCMD(cmd)

    # read in read hits
    R1Hits = []
    inFile = open(myData['R1Map'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)   
        R1Hits.append([pafLine['tStart'],pafLine['tEnd'],pafLine['strand']])        
    inFile.close()
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('R1Hits hits:\n')
        for v in R1Hits:
            fstream.write('%s\t%s\t%s\n' % (v[0],v[1],v[2]))
        fstream.flush()
    
    R2Hits = []
    inFile = open(myData['R2Map'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)   
        R2Hits.append([pafLine['tStart'],pafLine['tEnd'],pafLine['strand']])        
    inFile.close()
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('R2Hits hits:\n')
        for v in R2Hits:
            fstream.write('%s\t%s\t%s\n' % (v[0],v[1],v[2]))
        fstream.flush()

# now, ready to make plot.  we will use background of coverage info to help orient..

    myData['coverageInWindowsShowVectorESP'] = myData['bedWindowsFile'] + '.coverage.show-vector-esp.png'
    
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

    maxDepth = 1.1 * max(depths) # 10% more


    plt.figure(figsize=(8,6))        
    plt.plot(xpos,depths,'o',color='black',markersize=3)
    plt.xlim([1,myData['contigLen']]) 
    plt.ylim([0,maxDepth])
       
    plt.xlabel(myData['name'])
    plt.ylabel('Read Depth')   
    plt.title('Coverage with Vector and End Sequence Hits') 
    
    # show vector hits at depth 2, black for fwd, grey for reverse    
    for hit in vectorHits:
        b = hit[0] + 1
        e = hit[1]
        strand = hit[2]
        if strand == '+':
            myCol = 'black'
        else:
            myCol = 'grey'
            
        y = 2
        plt.plot([b,e],[y,y],color=myCol,linewidth=4)

    # show R1 hits   
    for hit in R1Hits:
        b = hit[0] + 1
        e = hit[1]
        strand = hit[2]
        if strand == '+':
            myCol = 'blue'
        else:
            myCol = 'cyan'            
        y = 2
        plt.plot([b,e],[y,y],color=myCol,linewidth=4)

    # show R1 hits   
    for hit in R2Hits:
        b = hit[0] + 1
        e = hit[1]
        strand = hit[2]
        if strand == '+':
            myCol = 'red'
        else:
            myCol = 'pink'            
        y = 2
        plt.plot([b,e],[y,y],color=myCol,linewidth=4)


    plt.savefig(myData['coverageInWindowsShowVectorESP'])
    plt.close()
#######################################################################    
def extract_segment(myData):
    faSeqs = read_fasta_file_to_dict(myData['contig'])
    
    if len(faSeqs) != 1:
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nERROR\nHave to many seqs in %s\n' % myData['contig'])
            fstream.flush()
        sys.exit()
    
    
    k = faSeqs.keys()
    k = list(k)
    seqName = k[0]
    seq = faSeqs[seqName]
    
    newSeq = seq['seq'][myData['trimStart']-1:myData['trimEnd']]    
    newSeqBreaks = add_breaks_to_line(newSeq)
    myData['newSeqFileName'] = myData['outDir'] + myData['name'] + '.fa'
        
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nExtracted length is:\t%i' % len(newSeq))
        fstream.write('\nwriting to new file:\t%s\n' % myData['newSeqFileName'])        
        fstream.flush()        
    outFile = open(myData['newSeqFileName'],'w')
    outFile.write('>%s\n' % myData['name'])
    outFile.write('%s\n' % (newSeqBreaks))
    outFile.close()
#######################################################################    
def run_racon(myData):
    myData['PAFout'] = myData['outDir'] + 'align.paf'
    myData['polishedFa'] = myData['outDir'] + myData['name'] + '.fa'

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
        fstream.write('\nminimap2 cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)        

    # racon
    cmd = 'racon --no-trimming %s %s %s > %s ' % (myData['longread'],myData['PAFout'],myData['contig'],myData['polishedFa'] )
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nracon cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)        
#######################################################################    
def run_rotate_and_remove(myData):
    
    myData['rotatedFa'] = myData['outDir'] + myData['name'] + '.rotated.fa'
    myData['rotatedCleanedFa'] = myData['outDir'] + myData['name'] + '.rotated.clean.fa'
    myData['finalFa'] = myData['outDir'] + myData['cloneName'] + '.fa'

    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nNOTE: Rotate based on mapping of vector.\n I am assuming there is a single circular contig\nand the the vector is not split\n')
        fstream.flush()
    
    
    myData['vectorMapOut'] = myData['outDir'] + 'vector-map.paf'
    
    cmd = 'minimap2 -c -x asm5 %s %s > %s ' % (myData['contig'],myData['vector'],myData['vectorMapOut'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nmap vector cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)

    # read in vector hits
    vectorHits = []
    inFile = open(myData['vectorMapOut'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)   
        vectorHits.append([pafLine['tStart'],pafLine['tEnd'],pafLine['strand'],pafLine['qStart'],pafLine['qEnd'],pafLine['qLen']])        
    inFile.close()
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('Vector hits:\n')
        for v in vectorHits:
            fstream.write('%s\t%s\t%s\t%s\t%s\n' % (v[0],v[1],v[2],v[3],v[4]))
        fstream.flush()
        
    if len(vectorHits) != 1:
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('ERROR! not 1 vector hit:\n')
            fstream.flush()
        sys.exit()    
    vectorHit = vectorHits[0]
    
    
    if vectorHit[2] == '-':
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('Need to make reverse complement of sequence!\n')
            fstream.flush()
    
        fastaSeqs = read_fasta_file_to_dict(myData['contig'])
        name = list(fastaSeqs.keys())[0]
        seq = fastaSeqs[name]['seq']
        seq = revcomp(seq)
        seq = add_breaks_to_line(seq,n=100)
        myData['assemFa'] = myData['outDir'] + myData['name']  + '.working.rc.fa'        
        outFile = open(myData['assemFa'],'w')
        outFile.write('>%s\n' % (name))
        outFile.write(seq)
        outFile.close()
    
        # then, need to redo the alignment
        myData['vectorMapOut'] = myData['outDir'] + 'vector-map.redo.paf'    
        cmd = 'minimap2 -c -x asm5 %s %s > %s ' % (myData['assemFa'],myData['vector'],myData['vectorMapOut'])
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nmap vector cmd is:\n%s\n' % cmd)
            fstream.flush()
        runCMD(cmd)
        # read in vector hits
        vectorHits = []
        inFile = open(myData['vectorMapOut'],'r')
        for line in inFile:
            line = line.rstrip()
            line = line.split()
            pafLine = parse_paf_line(line)   
            vectorHits.append([pafLine['tStart'],pafLine['tEnd'],pafLine['strand'],pafLine['qStart'],pafLine['qEnd'],pafLine['qLen']])        
        inFile.close()
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('Vector hits:\n')
            for v in vectorHits:
                fstream.write('%s\t%s\t%s\t%s\t%s\n' % (v[0],v[1],v[2],v[3],v[4]))
            fstream.flush()
        vectorHit = vectorHits[0]
    else:
        myData['assemFa'] = myData['contig']


    # check the size hit
    if vectorHit[3] != 0 or vectorHit[4] != vectorHit[5]:
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('ERROR! vector not totally spanned in hit:\n')
            fstream.flush()
        sys.exit()
    

    # for now, will assume that the vector sequence is not across the circle junction -- that would be a complication
    # will do the coordinates in 1 base system, not bed/psl
    tStart = vectorHit[0] + 1
    fastaSeqs = read_fasta_file_to_dict(myData['assemFa'])
    name = list(fastaSeqs.keys())[0]

    originalSeq = fastaSeqs[name]['seq']
    part1 = originalSeq[tStart-1:]  # starts at the original seq
    part2 = originalSeq[0:tStart-1]
    print('p1 len',len(part1))
    print('p2 len',len(part2))
    print('total',len(part1) + len(part2))

    outFile = open(myData['rotatedFa'],'w')
    outFile.write('>%s\n' % (myData['cloneName']))
    seq = part1 + part2
    seq = add_breaks_to_line(seq,n=100)
    outFile.write(seq)
    outFile.close()

    # now determine where vector is in the rotated
    myData['rotatedVectorMapOut'] = myData['outDir'] + 'rotated-vector-map.paf'
    cmd = 'minimap2 -c -x asm5 %s %s > %s ' % (myData['rotatedFa'],myData['vector'],myData['rotatedVectorMapOut'])
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nrotated map vector cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)

    # read in vector hits
    vectorHits = []
    inFile = open(myData['rotatedVectorMapOut'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)   
        vectorHits.append([pafLine['tStart'],pafLine['tEnd'],pafLine['strand'],pafLine['qStart'],pafLine['qEnd'],pafLine['qLen']])        
    inFile.close()
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('Vector hits:\n')
        for v in vectorHits:
            fstream.write('%s\t%s\t%s\t%s\t%s\n' % (v[0],v[1],v[2],v[3],v[4]))
        fstream.flush()
    vectorHit = vectorHits[0]
    
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nRemoving vector!:\n')
    fstream.flush()

    
    fastaSeqs = read_fasta_file_to_dict(myData['rotatedFa'])    
    name = list(fastaSeqs.keys())[0]
    seq = fastaSeqs[name]['seq']    
    #tEnd is in bedFormat -- so can use directly, will get the next base...
    newSeq = seq[vectorHit[1]:]

    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nvector seq ended at map pos %i\n' % vectorHit[1])
        fstream.write('length of new trimmed seq is %i\n' % len(newSeq))
        fstream.flush()
    
    newSeq = add_breaks_to_line(newSeq,n=100)
    outFile = open(myData['rotatedCleanedFa'],'w')
    outFile.write('>%s\n' % (myData['cloneName']))
    outFile.write(newSeq)
    outFile.close()
    

    cmd = 'cp %s %s' % ( myData['rotatedCleanedFa'],myData['finalFa'])    
    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\ncopy to final file cmd is:\n%s\n' % cmd)
        fstream.flush()
    runCMD(cmd)
    
#####################################################################
