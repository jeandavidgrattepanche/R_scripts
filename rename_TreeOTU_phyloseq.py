import sys
import os
import re
import time
import string
import os.path
from sys import argv
import numpy

#occurence = 13 (only number to sum after)

def Merged_taxa(treefile):
	out = open(treefile.split('.tree')[0]+"_phyloseq.tree","w+")
	with open(treefile) as f:
		for line in f:
			for elt in line.split('OTU'):
				if elt.split('_')[0].isnumeric() and elt.split('_')[1].split('r:')[0].isnumeric():
					rename = "OTU" + elt.split('_')[0]+':'+ elt.split('_')[1].split('r:')[1]
				else:
					rename = elt

				out.write(rename)

	out.close()
def main():
	script, treefile = argv
	Merged_taxa(treefile)
main()