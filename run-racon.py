import sys
import bactools3
import argparse
import os


parser = argparse.ArgumentParser(description='extract segment from contig')


parser.add_argument('--outdir', type=str,help='Path to output directory for this analysis',required=True)
parser.add_argument('--name', type=str,help='name of project',required=True)
parser.add_argument('--contig', type=str,help='fasta file of contig, should be single sequence',required=True)
parser.add_argument('--longread', type=str,help='file of long reads',required=True)



args = parser.parse_args()

# setup dictionary for holding information and passing to functions
myData = {}
myData['name'] = args.name 
myData['contig'] = args.contig 
myData['outDir'] = args.outdir 
myData['longread'] = args.longread 
myData['longreadtype'] = 'ont' # for now, only option is oxford nanopore



# setup needed files
if myData['outDir'][-1] != '/':
    myData['outDir'] += '/'

if os.path.isdir(myData['outDir']) is False:
    print('Output dir doest not exist, making it!')
    cmd = 'mkdir ' + myData['outDir']
    print(cmd, flush=True)
    bactools3.runCMD(cmd)

myData['logFileName'] = myData['outDir'] + 'run-racon.log'
print('logFileName',myData['logFileName'])


logi = 1
while os.path.isfile(myData['logFileName']) and os.path.getsize(myData['logFileName']) > 0:
    myData['logFileName'] = myData['outDir'] + 'run-racon.log.' + str(logi)
    logi += 1
    print('logFileName',myData['logFileName'])
    
myData['logFile'] = open(myData['logFileName'],'w')

bactools3.check_prog_paths(myData)
bactools3.write_initial_log(myData)

###########################################################################################
bactools3.run_racon(myData)




myData['logFile'].close()