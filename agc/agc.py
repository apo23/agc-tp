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
        

def get_chunks(sequence, chunk_size):
    chunk_num = int(len(sequence) / chunk_size)
    chunk_list = []
    if chunk_num < 4:
        chunk_num = int(len(sequence) / 4)
    count = 0
    for i in range(chunk_num):
        chunk_list.append(sequence[count:count+chunk_size])
        count += chunk_size
    return chunk_list


def test_unique():
    pass


def test_common():
    pass


def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence)-(kmer_size-1)):
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    kmers = cut_kmer(sequence, kmer_size)
    for kmer in kmers:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [id_seq]
        else:
            kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
    diffs = [x for x in alignment_list[0] if x != y for y in alignment_list[1]]
    ids = len(diffs)
    return ids/len(alignment_list[0])


def detect_chimera(perc_identity_matrix):
    stdevs = []
    diff_count = 0
    favored = -1
    prev_fav = -1
    switch_flag = 0
    for i in range(len(perc_identity_matrix)):
        curr_line = perc_identity_matrix[i]
        if prev_fav != favored:
            switch_flag += 1
        prev_fav = favored
        favored = [i for i, x in enumerate(curr_line) if x == max(curr_line)][0]
        curr_stdev = statistics.stdev(curr_line)
        stdevs.append(curr_stdev)
    if prev_fav != favored:
        switch_flag += 1
    if switch_flag >= 2 and statistics.mean(stdevs) > 5:
        return True
    return False


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    derep = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    seqs = []
    kmer_dict = {}
    id_seq = 0
    for seq, occ in derep:
        is_chimera = False
        tmp_id_mat = []
        tmp_chunks = get_chunks(seq, chunk_size)
        mates = search_mates(kmer_dict, sequence, kmer_size)
        if len(mates) >= 2:
            parent1 = seqs[mates[0]]
            parent2 = seqs[mates[1]]

            p1chunks = get_chunk(parent1[0], chunk_size)
            p2chunks = get_chunk(parent2[0], chunk_size)
            for c in range(len(tmp_chunks)): 
                line = []
                line.append(get_identity(global_align(seq, parent1[0])))
                line.append(get_identity(global_align(seq, parent2[0])))
                tmp_id_mat.append(line)
            is_chimera = detect_chimera(tmp_id_mat)
        kmer_dict = get_unique_kmer(kmer_dict, seq, id_seq, kmer_size
        if is_chimera == False:
            seqs.append([seq,occ])
            yield([seq, occ])
        


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
    test = dereplication_fulllength(args.amplicon_file, args.minseqlen, args.mincount)
    seqs = []
    count=0
    for t in test:
        print(count)
        count+=1
        seqs.append(t)
    print(len(seqs))
    
    


if __name__ == '__main__':
    main()
