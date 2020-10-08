import sys
import bactools3
import argparse
import os


print('Checking paths!\n')

myData = {}
bactools3.check_prog_paths(myData)

print('\nLooks ok, ready to run!')