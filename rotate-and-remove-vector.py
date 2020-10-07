import sys
import bactools3
import argparse
import os


parser = argparse.ArgumentParser(description='extract segment from contig')


parser.add_argument('--outdir', type=str,help='Path to output directory for this analysis',required=True)
parser.add_argument('--name', type=str,help='name of project',required=True)
parser.add_argument('--contig', type=str,help='fasta file of contig, should be single sequence',required=True)
parser.add_argument('--vector', type=str,help='fasta file of vector should be single sequence',required=True)
parser.add_argument('--clonename', type=str,help='clone name used to extract end sequences',required=True)


args = parser.parse_args()

# setup dictionary for holding information and passing to functions
myData = {}
myData['name'] = args.name 
myData['contig'] = args.contig 
myData['outDir'] = args.outdir 
myData['vector'] = args.vector
myData['cloneName'] = args.clonename 



# setup needed files
if myData['outDir'][-1] != '/':
    myData['outDir'] += '/'

if os.path.isdir(myData['outDir']) is False:
    print('Output dir doest not exist, making it!')
    cmd = 'mkdir ' + myData['outDir']
    print(cmd, flush=True)
    bactools3.runCMD(cmd)

myData['logFileName'] = myData['outDir'] + 'rotate-and-remove-vector.log'
print('logFileName',myData['logFileName'])


logi = 1
while os.path.isfile(myData['logFileName']) and os.path.getsize(myData['logFileName']) > 0:
    myData['logFileName'] = myData['outDir'] + 'rotate-and-remove-vector.log.' + str(logi)
    logi += 1
    print('logFileName',myData['logFileName'])
    
myData['logFile'] = open(myData['logFileName'],'w')

bactools3.check_prog_paths(myData)
bactools3.write_initial_log(myData)

###########################################################################################
bactools3.run_rotate_and_remove(myData)


myData['logFile'].close()