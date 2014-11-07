#! /usr/bin/python
# -*- coding: utf-8 -*-
#
# TODO:
# This script will work, if there are scheme.mres and Gene.gene files.

# cd folder/
# python TmiRUSite.py -s scheme.mres -g genedirectory -o resultname

from __future__ import division, print_function
import os
import sys
import re
import argparse


def arg_receive():
    "To take and get out the arguments"
    parser = argparse.ArgumentParser(
        description = 'The script allows for retrieving a set of \
        additional sequences at left and at right from the miRNA \
        binding site. This script search for binding sites of UTRs in \
        scheme results (scheme.mres) of the miRTarget program \
        and allows to search for mRNA sequence with miRNA binding sites. \
        The oligonucleotide contain 15 nucleotides before binding sites  \
        and 42 nucleotides after the start position. The GENE.gene \
        files used at searching for miRNA binding site.\
        \
        #Example of “res.mres” file scheme results  that has to folowing data:\
        #gene name	miRNA	number position	mRNA part	∆G	∆G/∆Gm	chromosome\
        #>SSTR3	>miR-147a	1	1085	CDS	-87.0	81.9	20\
        #5 - GCAUGAGCACCUGCCACAUGC - 3\
        #     ||| |||||||| ||||||||\
        #3 - CGUCUUCGUAAA-GGUGUGUG - 5\
        \
        #Example of “SSTR3.gene” file that has to folowing data:\
        #>SSTR3 | 1-525 | 526-1782 | 1783-2124\
        #CGCATCTCTCATCACTCCCCCTCATTCTGCCTTTCCTCCTACTCACGGTCTCCTCTC...',
        epilog = 'Send the bugs to devolia18@mail.ru',
        argument_default=argparse.SUPPRESS)
    
    parser.add_argument(
        '-s',
        type=str,
        required=True,
        nargs='?',
        metavar='The file name with schemes of binding sites scheme.mres',
        help='The file name with schemes of binding sites scheme.mres')
    parser.add_argument(
        '-g',
        type=str,
        required=True,
        nargs='?',
        metavar='The directory with GENE.gene files',
        help='The name of folder (directory) with GENE.gene files')
    parser.add_argument(
        '-o',
        type=str,
        required=True,
        nargs='?',
        metavar='Full text file with results',
        help='The file name with received results')
    
    args = parser.parse_args()
    
    return args.s, args.g, args.o


def read_site(scheme):
    dicgen = {}
    y = open(scheme, 'r')
    while True:
        x = y.readline()
        if not x:
            y.close()
            break
        else:
            if re.search('^>\w+', x):
                x = (x.rstrip()).split('\t')
                if 'CDS' not in x:
                    if '-v-' in x[0]:
                        xs = '-v-'
                    else:
                        xs = '-e-'
                    gname = ((x[0][1:]).split(xs))[0]
                    gpos = int(x[3])
                    gmir = x[1][1:]
                    domain = x[4]
                    inf = x[5:]
                    for i in inf:
                        inf[inf.index(i)] = float(i)
                    gstr = (y.readline()).rstrip()
                    gstr = (gstr.split(' '))[2]
                    gstr = gstr.lower()
                    if gname not in dicgen:
                        dicgen[gname] = {}
                    if gmir not in dicgen[gname]:
                        dicgen[gname][gmir] = [[gpos, inf, gstr, domain]]
                    else:
                        dicgen[gname][gmir].append([gpos, inf, gstr, domain])
    return(dicgen)

def main():
    
    scheme, gene_dir, res_file_name = arg_receive()
    
    scheme = os.path.abspath(scheme)
    genes = os.listdir(gene_dir)
    res_file_name = os.path.abspath(res_file_name)
    genes = filter(lambda x: x.endswith('.gene'), genes)
    dicgen = read_site(scheme)

    output = open(res_file_name, 'w')
    output.write('Target gene' + '\t' + 'miRNA' + '\t' + 'Position' + '\t' + 'Domain'+ '\t' + 'Energy' + '\t' + 'Score' + '\t' + 'Length' + '\t' + 'Binding site' + '\t' + 'mRNA fragment' + '\t' + '\n')
    for x in dicgen.keys():
        x2 = os.path.dirname(gene_dir) + '/' + x + '.gene'
        f2 = open(x2, 'r')
        g2 = (f2.readline()).split('|')
        if isinstance(g2, list):
            if '-' in g2[1]:
                g2[1].replace(' ', '')
                g2poscds = g2[1].split('-')
                g2poscds = int(g2poscds[0])
                if g2poscds == 0:
                    g2poscds += 1
            else:
                g2poscds = 1
        g2seq = ((f2.readline()).rstrip()).lower()
        g2seq = g2seq.replace('t', 'u')   
        f2.close()
        if isinstance(g2poscds, int):
            for y in dicgen[x]:
                for z in dicgen[x][y]:
                    pos = int(z[0])
                    pos = pos - g2poscds
                    energy = z[1]
                    bindseq = z[2]
                    domain = z[3]
                    amires = []
                    sites = []
                    if pos % 3 == 0:
                        poss = pos - 15
                        while poss <= (pos + 42):
                            k = g2seq[poss:(poss+3)] 
                            sites.append(k)
                            bs =  ''.join(sites)
                            amires.append(k)
                            poss += 3
                            l = ''.join(amires)
                        result_list = (x, y, pos, domain, energy, bindseq, bs)
                    elif ((pos - 1) % 3) == 0:
                        pos2v = pos - 16
                        while pos2v <= (pos + 41):
                            h = g2seq[pos2v:(pos2v+3)] 
                            sites.append(h)
                            bs =  ''.join(sites)
                            amires.append(h)
                            pos2v += 3
                            l = ''.join(amires)
                        result_list = (x, y, pos, domain, energy, bindseq, bs)

                    else:
                        pos3v = pos -17
                        while pos3v <= (pos + 40):
                            j = g2seq[pos3v:(pos3v+3)] 
                            sites.append(j)
                            bs =  ''.join(sites)
                            amires.append(j) 
                            pos3v += 3
                            l = ''.join(amires)
                        result_list = (x, y, pos, domain, energy, bindseq, bs)

                    for r in result_list:
                        if isinstance(r, list):
                            for r2 in r:
                                if isinstance(r2, list):
                                    for r3 in r:
                                        output.write(str(r3) + '\t')
                                else:
                                    output.write(str(r2) + '\t')
                        else:
                            output.write(str(r) + '\t')
                    output.write('\n')
    output.close()

main()
