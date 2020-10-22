#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
#import nwalign3 as nw

__author__ = "Apollinaire Roubert"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Apollinaire Roubert"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Apollinaire Roubert"
__email__ = "apollinaire.roubert@etu.u-paris,fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    with open(amplicon_file, "r") as filin:
        lines = filin.readlines()
        flag = 0
        tmp = ''
        for line in lines:
            if line.startswith(">"):
                if tmp != '':
                    if len(tmp) >= minseqlen:
                        yield tmp
                    tmp = ''
            else:
                flag = 0
                tmp += line.strip()
        if len(tmp) >= minseqlen:
            yield(tmp)
            


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dict_seq = read_fasta(amplicon_file, minseqlen)
    seqs = {}
    yield_list = []
    yield_index = []
    for seq in dict_seq:
        if seq not in seqs:
            seqs[seq] = 1
        else:
            seqs[seq] += 1
    for seq, occ in seqs.items():
        if occ >= mincount:
            yield_list.append(seq)
            yield_index.append(occ)
    while len(yield_index) != 0:
        curr = [i for i, x in enumerate(yield_index) if x == max(yield_index)]
        yield [yield_list[curr[0]], yield_index[curr[0]]]
        yield_list.pop(curr[0])
        yield_index.pop(curr[0])
        

def test_get_chunks():
    pass


def test_unique():
    pass


def test_common():
    pass


def test_cut_kmer():
    pass


def test_get_unique_kmer():
    pass


def test_search_mates():
    pass


def test_get_identity():
    pass


def test_detect_chimera():
    pass


def test_chimera_removal():
    pass


def test_abundance_greedy_clustering():
    pass


def test_write_OTU():
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    dereplication = dereplication_fulllength(args.amplicon_file, args.minseqlen,
                                    args.mincount)
    
    


if __name__ == '__main__':
    main()
