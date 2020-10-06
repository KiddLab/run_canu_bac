import sys
import bactools3
import argparse
import os


parser = argparse.ArgumentParser(description='run Canu on oxford nanopore BACs')


parser.add_argument('--outDirBase', type=str,help='Path to output directory',required=True)
parser.add_argument('--name', type=str,help='name of project',required=True)
parser.add_argument('--longread', type=str,help='file of long reads',required=True)
parser.add_argument('--contam', type=str,help='fasta of contamination (E. coli), with bwa mem index',required=True)
parser.add_argument('--target',type=str,help='fasta file of sequnece to use for begining of rotated sequence',required=True)
parser.add_argument('--threads',type=int,default=3,help='Number of threads to use, default is 3')


args = parser.parse_args()

# setup dictionary for holding information and passing to functions
myData = {}
myData['outDirBase'] = args.outDirBase 
myData['name'] = args.name 
myData['longread'] = args.longread 
myData['longreadtype'] = 'ont' # for now, only option is oxford nanopore
myData['contam'] = args.contam 
myData['numThreads'] = args.threads 
myData['targetFa']  = args.target 

# setup needed files
if myData['outDirBase'][-1] != '/':
    myData['outDirBase'] += '/'

if os.path.isdir(myData['outDirBase']) is False:
    print('Output dir doest not exist, making it!')
    cmd = 'mkdir ' + myData['outDirBase']
    print(cmd, flush=True)
    bactools3.runCMD(cmd)

myData['logFileName'] = myData['outDirBase'] + 'run-initial-assembly.log'
print('logFileName',myData['logFileName'])


logi = 1
while os.path.isfile(myData['logFileName']) and os.path.getsize(myData['logFileName']) > 0:
    myData['logFileName'] = myData['outDirBase'] + 'run-initial-assembly.log.' + str(logi)
    logi += 1
    print('logFileName',myData['logFileName'])
    
myData['logFile'] = open(myData['logFileName'],'w')



print('Will run with %i threads!' % myData['numThreads'], flush=True)
bactools3.check_prog_paths(myData)
bactools3.write_initial_log(myData)

# step 1 map long read data against e coli and filter out hits
bactools3.filter_contam_longread(myData)
bactools3.print_longread_stats(myData['longread'],myData)
bactools3.print_longread_stats(myData['longReadFilt'],myData)

# step 2 run canu assembly
bactools3.run_canu_assem(myData)


###############################################################################

# cleanup
myData['logFile'].write('\nEXIT\n')
now = datetime.datetime.now()
myData['logFile'].write(str(now))
myData['logFile'].write('\n')
myData['logFile'].close()


