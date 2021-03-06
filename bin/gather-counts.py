#! /usr/bin/env python
"""
This script gathers & converts Salmon output counts into something that
edgeR can read ("counts files").

Run it in a directory above all of your Salmon output directories, and
it will create a bunch of '.counts' files that you can load into R.

See https://github.com/ngs-docs/2015-nov-adv-rna/ for background info.

C. Titus Brown, 11/2015
"""
import os, os.path
import sys
import csv

def process_quant_file(root, filename, outname, dirname):
    """
    Convert individual quant.sf files into .counts files (transcripts\tcount).
    """
    # print >>sys.stderr, 'Loading counts from:', root, filename
    outfp = open(os.path.join(root, outname), 'w')
    outfp.write("transcript\t{}\n".format(dirname))

    d = {}
    full_file = os.path.join(root, filename)
    for line in open(full_file):
        if line.startswith('Name'):
            continue
        name, length, eff_length, tpm, count = line.strip().split('\t')

        outfp.write("{}\t{}\n".format(name, float(tpm)))


def main():
    """
    Find all the quant.sf files, convert them into properly named .counts
    files.

    Here, "proper name" means "directory.counts".
    """
    quantlist = []

    start_dir = sys.argv[1]
    for root, dirs, files in os.walk(start_dir):
        for filename in files:
            if filename.endswith('quant.sf'):
                dirname = os.path.basename(root)
                outname = dirname + '.tpm'
                process_quant_file(root, filename, outname, dirname)
                quantlist.append(outname)
                
                break

    #print ",\n".join([ "\"%s\"" % i for i in sorted(quantlist)])

if __name__ == '__main__':
    main()
