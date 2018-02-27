#!/usr/bin/env python3

'''
Below are the global variables and dicts for Clinvar. Modify variables
with caution, as the program can become unstable.

Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

__version__ = "0.2.3"


### getncbi.py global variables 

# batch size for querying NCBI, not above 5000
batch_size = 4500

# name of software for ncbi query tagging (preferred by NCBI)
etool = 'Clinotator'


### vcf.py global variables 

# limit for the vaf header lines to search. They are usually <100, 200 is safe
max_vcf_header_size = 200

# header INFO fields to be inserted into the output vcf
new_headers = ['##INFO=<ID=VID,Number=1,Type=Integer,Description="ClinVar var'
               'iation ID">\n',
               '##INFO=<ID=CVVT,Number=A,Type=String,Description="ClinVar var'
               'iant type">\n',
               '##INFO=<ID=CVMA,Number=A,Type=String,Description="ClinVar min'
               'or allele">\n',
               '##INFO=<ID=CVCS,Number=A,Type=String,Description="ClinVar cli'
               'nical significance">\n',
               '##INFO=<ID=CVSZ,Number=A,Type=Integer,Description="ClinVar st'
               'ars">\n',
               '##INFO=<ID=CVNA,Number=A,Type=Integer,Description="ClinVar nu'
               'mber of clinical assertions">\n',
               '##INFO=<ID=CVDS,Number=A,Type=String,Description="ClinVar con'
               'ditions">\n',
               '##INFO=<ID=CVLE,Number=A,Type=String,Description="ClinVar las'
               't evaluated">\n',
               '##INFO=<ID=CTRS,Number=A,Type=Float,Description="Clinotator r'
               'aw score">\n',
               '##INFO=<ID=CTAA,Number=A,Type=Float,Description="Clinotator a'
               'verage clinical assertion age">\n',
               '##INFO=<ID=CTWS,Number=A,Type=String,Description="Clinotator '
               'weighted significance">\n',
               '##INFO=<ID=CTRR,Number=A,Type=String,Description="Clinotator '
               'reclassification recommendation">\n']


### variation.py global vars

# dict of the different ClinVar variant star classifications
star_dict = {'practice guideline': 4,
             'reviewed by expert panel': 3,
             'criteria provided, multiple submitters, no conflicts': 2,
             'criteria provided, conflicting interpretations': 1,
             'criteria provided, single submitter': 1,
             'no assertion for the individual variant': 0,
             'no assertion criteria provided': 0,
             'no assertion provided': 0}        

# dict of reviewer status weights for each assertion
cutoff = {'practice guideline': 1.25,
          'reviewed by expert panel': 1.10,
          'criteria provided, single submitter': 1.0,
          'no assertion for the individual variant': 0.0,
          'no assertion criteria provided': 0.0,
          'no assertion provided': 0.0}

# dict of assertion significances for scoring
significance = {'Benign': -5,
                'Likely benign': -3,
                'Uncertain significance': -0.3,
                'Likely pathogenic': 1.6,
                'Pathogenic': 2.9,
                'drug response': 0, 'association': 0, 'risk factor': 0,
                'protective': 0, 'Affects': 0,
                'conflicting data from submitters': 0, 'other': 0,
                'not provided': 0}

# list of weighted score upper bounds for ctws bins
ctws_cutoffs = [('Benign', -25.8),
                ('Benign/Likely benign', -6),
                ('Likely benign', -3.84),
                ('Uncertain significance', 1.01),
                ('Likely pathogenic', 3.2),
                ('Pathogenic/Likely pathogenic', 14.17),
                ('Pathogenic', 10000000)]
