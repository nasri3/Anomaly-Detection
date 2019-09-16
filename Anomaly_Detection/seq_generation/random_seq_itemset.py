#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Random sequence generator.

It includes simple random generators of sequence of items or itemsets (not based on patterns).

"""

__author__ = "Thomas Guyet"
__copyright__ = "Copyright 2019, AGROCAMPUS-OUEST/IRISA"
__license__ = "LGPL"
__version__ = "1.0.1"
__maintainer__ = "Thomas Guyet"
__email__ = "thomas.guyet@irisa.fr"
    

import seqdb_generator as seqgen
import numpy as np


def gen_seqdb(seqnb, length, maxitems, fixedlength):
    """Fully random sequence generation. Generate sequences of itemsets.

    The function generates the file.

    :param seqnb: number of sequences
    :param length: mean length of the sequences
    :param maxitems: size of the vocabulary (maximum number of items)
    :param outputfile: output filename
    :return: generated sequences
    """
    sequences=[]
    item_gen = seqgen.item_generator(maxitems)
    for i in range(seqnb):
        s= seqgen.sequence()
        if fixedlength:
            rlength=length
        else:
            rlength=int(np.random.normal(length,max(1,int(length/10))))
    
        for j in range(rlength):
            s.seq.append( (item_gen.generate(),j) )
    
        sequences.append(s)
    return sequences
    
def gen_seq_itemset():

    fout = open(outputfile, "w")
    for i in range(seqnb):
        if i!=0:
            fout.write("\n")
        if fixedlength:
            rlength=length
        else:
            rlength=int(np.random.normal(length,max(1,int(length/10))))
        for j in range(rlength):
            alreadyAdded=[]
            itemset=[]
            for k in range(itemsetSize):
                item = np.random.randint(0,maxitems+1,1)[0]
                while item in alreadyAdded:
                    item = np.random.randint(0,maxitems+1,1)[0]
                alreadyAdded.append(item)
                itemset.append(item)
            itemset.sort()
            for it in itemset[:-1]:
                if asp:
                    fout.write( "seq("+str(i) +","+str(j) +","+str(it) +"). ")
                else:
                    fout.write( str(it) + ":")
            if asp:
                fout.write( "seq("+str(i) +","+str(j) +","+str(itemset[len(itemset)-1]) +"). ")
            else:
                fout.write( str(itemset[len(itemset)-1]) )
            if j!=(rlength-1) and not asp:
                fout.write(",")
    fout.close()
