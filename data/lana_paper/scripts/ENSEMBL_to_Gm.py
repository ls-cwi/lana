#!/usr/bin/python
"""converts ENSEMBL orthologs to candidate list format

Author:     Gunnar W. Klau
Version:    0.9
"""

import sys
import re
import os
import glob
from optparse import OptionParser
import traceback
from collections import OrderedDict

parser = OptionParser()
parser.add_option("-i", dest="ENSEMBL_file", help="ENSEMBL input", metavar="FILE")
parser.add_option("-o", dest="cand_list_file", help="candidate list file", metavar="FILE")
parser.add_option("--s1", dest="species_string_1", help="optional species ID for first protein in lists", type="string")
parser.add_option("--s2", dest="species_string_2", help="optional species ID for remaining proteins in lists", type="string")

(options, args) = parser.parse_args()

if not options.ENSEMBL_file or not options.cand_list_file:   # if filename is not given
    parser.error('input or output file not given')

# we keep the orthologs in a dictionary
O = {}

i = open(options.ENSEMBL_file)
o = open(options.cand_list_file, 'w')

print "parsing", options.ENSEMBL_file, "...",
sys.stdout.flush()
line = i.readline() # eat header line
while line:
    line = i.readline()
    L = line.split()
    if len(L) == 4: # these are the lines with orthologs, L[1] and L[3] are the ENSEMBL protein IDs
      if not O.get(L[1]): O[L[1]] = []
      O[L[1]].append(L[3])
print "done."

print "writing", options.cand_list_file, "...",
for mp in O:
    if options.species_string_1: o.write(options.species_string_1)
    o.write(mp)
    O[mp] = list(OrderedDict.fromkeys(O[mp])) # make unique, gthere were some doubles in the list...
    for hp in O[mp]:
        o.write(' ')
        if options.species_string_2: o.write(options.species_string_2)
        o.write(hp)
    o.write("\n")
print "done."

#compute some stats
hist = {}
for mp in O:
    if not hist.get(len(O[mp])): hist[len(O[mp])] = 1
    else: hist[len(O[mp])] = hist[len(O[mp])] + 1

print "some stats (histogram): ", hist

i.close()
o.close()
