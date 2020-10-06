import sys
import bactools3
import argparse
import os


parser = argparse.ArgumentParser(description='make images of contig')


parser.add_argument('--outdir', type=str,help='Path to output directory for this analysis',required=True)
parser.add_argument('--name', type=str,help='name of project',required=True)
parser.add_argument('--longread', type=str,help='file of long reads',required=True)
parser.add_argument('--contig', type=str,help='fasta file of contig, should be single sequence',required=True)

args = parser.parse_args()

# setup dictionary for holding information and passing to functions
myData = {}
myData['name'] = args.name 
myData['longread'] = args.longread 
myData['longreadtype'] = 'ont' # for now, only option is oxford nanopore
myData['contig'] = args.contig 
myData['outDir'] = args.outdir 

# setup needed files
if myData['outDir'][-1] != '/':
    myData['outDir'] += '/'

if os.path.isdir(myData['outDir']) is False:
    print('Output dir doest not exist, making it!')
    cmd = 'mkdir ' + myData['outDir']
    print(cmd, flush=True)
    bactools3.runCMD(cmd)

myData['logFileName'] = myData['outDir'] + 'make-assem-images.log'
print('logFileName',myData['logFileName'])


logi = 1
while os.path.isfile(myData['logFileName']) and os.path.getsize(myData['logFileName']) > 0:
    myData['logFileName'] = myData['outDir'] + 'make-assem-images.log.' + str(logi)
    logi += 1
    print('logFileName',myData['logFileName'])
    
myData['logFile'] = open(myData['logFileName'],'w')

bactools3.write_initial_log(myData)

###########################################################################################
# run miropeats self
bactools3.run_miropeats_self(myData)

# get coverage and long-read information
# map to contig
bactools3.map_to_contig_paf(myData)

bactools3.make_windows_bed(myData,500,100)
bactools3.make_coverage_plot(myData)

bactools3.make_gc_plot(myData)
# show positions of 5 longest reads
bactools3.make_coverage_plot_showlong(myData,10)


myData['logFile'].close()