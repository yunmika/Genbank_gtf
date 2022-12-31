#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Description :   This is a python script about genbank to gtf
@Author      :   Fuchuan Han 
@Time        :   2022/09/12 23:38:29
"""


from operator import index
from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse
import csv
import os


def Argparse():
    '''
    Adding parameters
    '''
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-i', '--input', help='Please input the genbank file')
    parser.add_argument('-o', '--output', help='output file')
    return parser.parse_args()

def ToGtf():
    '''
    Read the file and generate the gtf format list
    '''
    gtf_lt = []
    # read the mitogenome genbank file
    for records in SeqIO.read(args.input, 'genbank').features:
        if records.type == 'CDS':
            # determine if only one CDS is included
            if records.location_operator is None:
                gtf_lt.append(['Chr1', 'Biopy', 'CDS',
                            # the start location of CDS
                            int(records.location.start) + 1, 
                            # the end location of CDS
                            records.location.end, 
                            '.', 
                            '+' if records.strand == 1 else '-', 
                            # the codon_start loaction
                            records.qualifiers['codon_start'][0],
                            'gene_id \"{0}\"; transcript_id \"{0}\";'.format(records.qualifiers['gene'][0])])
                # Copy the CDS line definition as an exon
                gtf_lt.append(['Chr1', 'Biopy', 'exon',
                            # the start location of CDS
                            int(records.location.start) + 1, 
                            # the end location of CDS
                            records.location.end, 
                            '.', 
                            '.', 
                            # the codon_start loaction
                            records.qualifiers['codon_start'][0],
                            'gene_id \"{0}\"; transcript_id \"{0}\";'.format(records.qualifiers['gene'][0])])
        elif records.type == 'exon':
            gtf_lt.append(['Chr1', 'Biopy', 'exon', 
                        int(records.location.start) + 1,
                        records.location.end,
                        '.',
                        '+' if records.strand == 1 else '-',
                        '.',
                        'gene_id \"{0}\"; transcript_id \"{0}\";'.format(records.qualifiers['gene'][0])])
    return gtf_lt

if __name__ == '__main__':
    args = Argparse()
    gtf_dat = ToGtf()
    gtf_ary = np.array(gtf_dat)
    gtf_dtf = pd.DataFrame(gtf_ary, index=None)
    # setting quoting=csv.QUOTE_NONE, quotechar=None (delete the "")
    gtf_dtf.to_csv('%s'%(args.output), sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, quotechar=None)
    os.system('sort -n -k 4 {0} > {1}'.format(args.output, 'sort_%s'%args.output))
