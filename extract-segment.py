import sys
import bactools3
import argparse
import os


parser = argparse.ArgumentParser(description='extract segment from contig')


parser.add_argument('--outdir', type=str,help='Path to output directory for this analysis',required=True)
parser.add_argument('--name', type=str,help='name of project',required=True)
parser.add_argument('--contig', type=str,help='fasta file of contig, should be single sequence',required=True)
parser.add_argument('--clonename', type=str,help='clone name used to extract end sequences',required=True)
parser.add_argument('--usecanuseg', action='store_true')



args = parser.parse_args()

# setup dictionary for holding information and passing to functions
myData = {}
myData['name'] = args.name 
myData['contig'] = args.contig 
myData['outDir'] = args.outdir 
myData['cloneName'] = args.clonename 
myData['useCanuSegment'] = args.usecanuseg 


# setup needed files
if myData['outDir'][-1] != '/':
    myData['outDir'] += '/'

if os.path.isdir(myData['outDir']) is False:
    print('Output dir doest not exist, making it!')
    cmd = 'mkdir ' + myData['outDir']
    print(cmd, flush=True)
    bactools3.runCMD(cmd)

myData['logFileName'] = myData['outDir'] + 'extract-segment.log'
print('logFileName',myData['logFileName'])


logi = 1
while os.path.isfile(myData['logFileName']) and os.path.getsize(myData['logFileName']) > 0:
    myData['logFileName'] = myData['outDir'] + 'extract-segment.log.' + str(logi)
    logi += 1
    print('logFileName',myData['logFileName'])
    
myData['logFile'] = open(myData['logFileName'],'w')

bactools3.write_initial_log(myData)

###########################################################################################
# get coordinates to extract
if myData['useCanuSegment'] is True:
    inFile = open(myData['contig'],'r')
    headerLine = inFile.readline()
    headerLine = headerLine.rstrip()
    inFile.close()    

    for fstream in [sys.stdout,myData['logFile']]:
        fstream.write('\nread in header line\n%s\n' % (headerLine) )
        fstream.flush()
        
    trimParts = headerLine.split()[-1]
    # check if trim is there at the end
    if 'trim=' != trimParts[0:5]:
        for fstream in [sys.stdout,myData['logFile']]:
            fstream.write('\nERRORRR! Did not find trim= at end')
            fstream.flush()
        sys.exit()
        
    p = trimParts.split('=')[1].split('-')
    b = int(p[0])
    e = int(p[1])
     
    myData['trimStart'] = b  # note, is 1 based
    myData['trimEnd'] = e  # note, is 1 based
else:
    sys.exit() # to implement
    
for fstream in [sys.stdout,myData['logFile']]:
    fstream.write('\nwill trim to %i - %i\n' % (myData['trimStart'],myData['trimEnd']) )	
    fstream.flush()
    


bactools3.extract_segment(myData)


            


    
    

    




###########################################################################################