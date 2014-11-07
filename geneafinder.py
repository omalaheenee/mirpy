#! /usr/bin/python
# -*- coding: utf-8 -*-

# To do (in terminal):
# cd folder/
# python GeneAFinder.py -g ~/.../gene_dictionary.txt/ -f ~/.../abstract.txt > file.txt
# ~/.../gene_dictionary.txt/      is directory to gene dictionary file (gene_dictionary.txt)
# ~/.../abstract.txt              is directory to file with abstracts (abstract.txt)  
# "File.txt" is file with results.
# The keywords can be changed in the script.

from __future__ import division, print_function
import os
import sys
import re
import argparse

def arg_receive():
    "This script To take and get out the arguments"
    parser = argparse.ArgumentParser(
        description = 'For use of this script it is necessary have \
        "abstract.txt" file with abstacts (text format of PubMed \
        (http://www.ncbi.nlm.nih.gov/pubmed/)). There are must be 2 \
        empty strings between abstracts. The result of searching for \
        contents of the following data: PMID (publication medicine \
        identification number) of the article,\ miRNA name, disease \
        localization (organ or tissue), keywords (methods, change fold,\
        cellular processes, functions, animal species, types of cells,\
        biological liquids, etc.) and genes. The results will be \
        selected according to each found gene in the abstract.',
        epilog = 'Send the bugs to devolia18@mail.ru',
        argument_default=argparse.SUPPRESS)
    
    parser.add_argument(
        '-g',
        type=str,
        required=True,
        nargs='?',
        metavar='The "gene_dictionary.txt" file',
        help='The name of file with genes ("gene_dictionary.txt")')

    parser.add_argument(
        '-f',
        type=str,
        required=True,
        nargs='?',
        metavar='The "abstract.txt" file',
        help='The name of file with abstracts ("abstract.txt")')
        
    args = parser.parse_args()
    
    return args.g, args.f


def article(txt_file_name, gene_dict):

    re_mir = re.compile(r'(miR|let)\-\d+\w*(-3p|-5p)?\*?', re.IGNORECASE)
    re_words = re.compile(r'\b([A-Z0-9]{1,})?(-)?([A-Z0-9]{1,})\b', re.VERBOSE)
    re_dwords = re.compile(r"""breast|
                            mammary|
                            lung|
                            (gastr|thym|pancreat)ic|
                            ((gastro)?intestin|gynecologic|nasopharynge|neuroectoderm|endometri|cervic|o?esophage|or|(colo)?rect)al|
                            endometriosis|
                            colon|
                            stomach|
                            liver|
                            hepat(ic|ocellular)|
                            prostate|
                            testi(s|(cular))|
                            ovar(ian|y)|
                            bowel|
                            intestine|
                            bladder|
                            kidney|
                            renal|
                            neck|
                            head|
                            brain|
                            bone|
                            thyroid|
                            skin|
                            leuk(a)?emia|
                            (bili|pituit)ary""", re.IGNORECASE | re.VERBOSE)
    re_kwords = re.compile(r"""immune|
                            blood|
                            serum|
                            saliva|
                            ((stem)(\ *))?cell(((\ *)|(-))(cycle|line))?|
                            tumo(u)?r(\ *)suppressor|
                            oncogene|                            
                            (metasta|apopto)sis|
                            (in|de)crease|
                            \b(up)|
                            \b(down)|
                            \b(low)|
                            \b(high)|
                            (\w+)(\ *)syndrome|
                            virus|
                            hypoxia|
                            neoplasm|
                            (T|B)-cell|
                            human|
                            hsa-|
                            mouse|
                            rat """, re.IGNORECASE | re.VERBOSE)
    re_swords = re.compile(r"""subtype|
                            basal(-like)?|
                            luminal((\ *)(a|b|c))?|
                            her2|
                            triple|
                            (N)?SCLC|
                            ((non-)?small|large)(\ *)cell|
                            squamous|
                            bronchioloalveolar|
                            stromal|
                            \b(ewing)\b|
                            \b(GIST)\b|
                            ileu(m|s)|
                            (cec|duoden|jejun)um|
                            appendix|
                            \w+oma(s)?\b""", re.IGNORECASE | re.VERBOSE)
    re_pat_s = re.compile(r'\s+')

    txt_file_name = os.path.abspath(txt_file_name)
    
    txt_file = open(txt_file_name, 'r')
    abstract = ''
    dejavu = []
    result_dict = {}
    while True:
        string = txt_file.readline()
        if not string:
            txt_file.close()
            break

        if not string.startswith('PMID'):
            abstract += string
        else:
            abstract_pmid = re_pat_s.split(string.rstrip())[1]
            result_dict[abstract_pmid] = {'word': [],
                                          'miR': [], 
                                          'dword': [],
                                          'sword': [],
                                          'kword': []}
            resw = re_words.finditer(abstract, re.MULTILINE)
            for r in resw:
                r = str(r.group(0))
                if r in gene_dict and r not in dejavu:
                    dejavu.append(r)
                    result_dict[abstract_pmid]['word'].append(r + "\t" + gene_dict[r])

            resm = re_mir.finditer(abstract, re.MULTILINE)
            for r in resm:
                r = str(r.group(0)).lower()
                if r not in dejavu:
                    dejavu.append(r)
                    result_dict[abstract_pmid]['miR'].append(r)

            resw = re_dwords.finditer(abstract, re.MULTILINE)
            for r in resw:
                r = str(r.group(0)).lower()
                if r not in dejavu:
                    dejavu.append(r)
                    result_dict[abstract_pmid]['dword'].append(r)

            resw = re_kwords.finditer(abstract, re.MULTILINE)
            for r in resw:
                r = str(r.group(0)).lower()
                if r not in dejavu:
                    dejavu.append(r)
                    result_dict[abstract_pmid]['kword'].append(r)

            resw = re_swords.finditer(abstract, re.MULTILINE)
            for r in resw:
                r = str(r.group(0)).lower()
                if r not in dejavu:
                    dejavu.append(r)
                    result_dict[abstract_pmid]['sword'].append(r)

            abstract = ''
            dejavu = []

    return result_dict

if __name__ == '__main__':
    
    gene_dict_file, article_file = arg_receive()
    gene_dict = {}
    for x in open(os.path.abspath(gene_dict_file)):
        x = x.rstrip().split('\t')
        gene_dict[x[0]] = ''
        for y in x[1:-1]:
            gene_dict[x[0]] += y
    
    print('PMID', 'Gene', 'Full gene name', 'Gene chromosome', 'Gene disease', 'miRNA', 'Cancer', 'Subtype', 'Feature', sep='\t')
    
    rd = article(article_file, gene_dict)
    for x1 in rd.keys():
        for x2 in rd[x1]['word'] :
                print (x1, x2, ', '.join(rd[x1]['miR']),
                               ', '.join(rd[x1]['dword']),
                               ', '.join(rd[x1]['sword']),
                               ', '.join(rd[x1]['kword']), sep='\t')
