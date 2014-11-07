#! /usr/bin/python
# -*- coding: utf-8 -*-

# TODO:
# This script will work, if there are scheme.mres and Gene.gene files.

# For install
# cd folder/
#  python TmiROSite.py -s scheme.mres -g directorygenes -o resultname

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
        binding site. This script search for binding sites of CDS in \
        scheme results (scheme.mres) of the miRTarget program \
        and allows to search for mRNA sequence with miRNA binding sites. \
        mRNA sequence with binding site and its encoded for oligopeptide \
        in open reading frame (ORF) can be defined by TmiROSite script.\
        The oligopeptides contain five amino acids before binding site \
        that is equivalent to 15 nucleotide and some amino acids after them.\
        It translate 42 nucleotides after start position of binding site. \
        The quantity of additional nucleotides can be increased if necessary.\
        The GENE.gene files used at searching for miRNA binding site. \
        Stop codons are designed by dash.\
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
        metavar='The file name with schemes of binding sites res.mres',
        help='The file name with schemes of binding sites res.mres')
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
                if 'CDS' in x:
                    if '-v-' in x[0]:
                        xs = '-v-'
                    else:
                        xs = '-e-'
                    gname = ((x[0][1:]).split(xs))[0]
                    gpos = int(x[3])
                    gmir = x[1][1:]
                    inf = x[5:]
                    gstr = (y.readline()).rstrip()
                    gstr = (gstr.split(' '))[2]
                    gstr = gstr.lower()
                    if gname not in dicgen:
                        dicgen[gname] = {}
                    if gmir not in dicgen[gname]:
                        dicgen[gname][gmir] = [[gpos, inf, gstr]]
                    else:
                        dicgen[gname][gmir].append([gpos, inf, gstr])
                    
    return(dicgen)

def main():
    amino_dict = {
        'uuu': 'F', 'uuc': 'F',
        'uua': 'L', 'uug': 'L', 'cuu': 'L', 'cuc': 'L', 'cua': 'L', 'cug': 'L',
        'ucu': 'S', 'ucc': 'S', 'uca': 'S', 'ucg': 'S',  'agu': 'S',  'agc': 'S',
        'uau': 'Y', 'uac': 'Y',
        'uaa': '-', 'uag': '-', 'uga': '-',
        'ugu': 'C', 'ugc': 'C',
        'ugg': 'W',
        'auu': 'I', 'auc': 'I', 'aua': 'I',
        'aug': 'M',
        'guu': 'V', 'guc': 'V', 'gua': 'V', 'gug': 'V',
        'ccu': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P',
        'acu': 'T', 'acc': 'T', 'aca': 'T', 'acg': 'T',
        'gcu': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A',
        'cau': 'H', 'cac': 'H',
        'caa': 'Q', 'cag': 'Q',
        'aau': 'N', 'aac': 'N',
        'aaa': 'K', 'aag': 'K',
        'gau': 'D', 'gac': 'D',
        'gaa': 'E', 'gag': 'E',
        'cgu': 'R', 'cgc': 'R', 'cga': 'R', 'cgg': 'R', 'aga': 'R', 'agg': 'R',
        'ggu': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G',
    }

    scheme, gene_dir, res_file_name = arg_receive()
    
    scheme = os.path.abspath(scheme)
    genes = os.listdir(gene_dir)
    res_file_name = os.path.abspath(res_file_name)
    genes = filter(lambda x: x.endswith('.gene'), genes)
    dicgen = read_site(scheme)


    output = open(res_file_name, 'w')
    output.write('Target gene' + '\t' + 'miRNA' + '\t' + 'Position' + '\t' + 'Energy' + '\t' + 'Score' + '\t' + 'Length' + '\t' + 'Binding site' + '\t' + 'Oligonucleotide' + '\t' + 'Oligopeptide' + '\t' + 'RF1' + '\t' + 'RF2' + '\t' + 'RF3' + '\n')
    for x in dicgen.keys():
        x2 = os.path.dirname(gene_dir) + '/' + x + '.gene'
        f2 = open(x2, 'r')
        g2 = (f2.readline()).split('|')
        if isinstance(g2, list):
            if '-' in g2[1]:
                g2[1].replace(' ', '')
                g2poscds = g2[1].split('-')
                g2poscds = int(g2poscds[-1])
        g2seq = ((f2.readline()).rstrip()).lower()
        g2seq = g2seq.replace('t', 'u')
        f2.close()

        for y in dicgen[x]:
            for z in dicgen[x][y]:
                pos = int(z[0]) - 1
                posnew = pos - g2poscds
                energy = z[1]
                ss = []
                bindseq = z[2]
                    
                for yy in range(0, 3):
                    zz = []
                    if yy in range(1,3):
                            zz.append("_")
                    while (yy+2) <= len(bindseq):
                        xx = bindseq[yy:(yy+3)]
                        yy = yy + 3
                        if xx in amino_dict:
                            zz.append(amino_dict[xx])
                        else:
                            zz.append("_")
                    ss.append(''.join(zz))
                    
                amires = []
                sites = []
                if posnew % 3 == 0:
                    pos1v = pos - 15
                    while pos1v <= (pos + 42):
                        k = g2seq[pos1v:(pos1v+3)]
                        sites.append(k)
                        if k in amino_dict:
                            amires.append(amino_dict[k])
                        pos1v += 3
                elif ((posnew - 1) % 3) == 0:
                    pos2v = pos - 16
                    while pos2v <= (pos + 41):
                        h = g2seq[pos2v:(pos2v+3)]
                        sites.append(h)
                        if h in amino_dict:
                            amires.append(amino_dict[h])
                        pos2v += 3
                elif ((posnew - 2) % 3) == 0:
                    pos3v = pos - 17
                    while pos3v <= (pos + 40):
                        j = g2seq[pos3v:(pos3v+3)]
                        sites.append(j)
                        if j in amino_dict:
                            amires.append(amino_dict[j])
                        pos3v += 3
                result_list = (x, y, pos+1, energy, bindseq, ''.join(sites), ''.join(amires), ss)
                

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
