#!/usr/bin/env python3

'''
Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
'''

import argparse, logging, sys
import pandas as pd
import bin.getncbi as getncbi
import bin.vcf as vcf
from Bio import Entrez

__version__ = "0.1.0"

# error logging function plus config
def error_handling():
    return ' {}. {}, line: {}'.format(sys.exc_info()[0],
                                           sys.exc_info()[1],
                                           sys.exc_info()[2].tb_lineno)

# argparse function
def getargs():
    # import argparse
    
    parser = argparse.ArgumentParser(prog='clinotator', formatter_class=argparse.RawTextHelpFormatter,
                                     description='Clinical interpretation of ambiguous ClinVar annotations')
    parser.add_argument("--version", action='version', version='\n'.join(['Clinotator v'+__version__, __doc__]))
    parser.add_argument("input", metavar=('file'), help="input file(s) (returns outfile for each)",
                        nargs='+')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-t', dest='type', choices=['vid', 'rsid', 'vcf'], required=True,
                        help='vid - ClinVar Variation ID list\n'
                        'rsid - dbSNP rsID list\n'
                        'vcf - vcf file (output vcf generated)')
    requiredNamed.add_argument('-e', dest='email', required=True,
                               help='NCBI requires an email for querying their databases')
    return parser.parse_args()

# # # # how to handle file types 
def input_selection(type, file, vid_list):
    try:
        if type == 'vid':
            with open(file) as f:
                vid_list = [line.rstrip('\n') for line in f]
            return vid_list
        
        elif type == 'rsid':
            with open(file) as f:
                rsid_list = [line.lstrip('rsRS').rstrip('\n') for line in f]
                getncbi.rsid_to_vid(rsid_list, vid_list)
            return vid_list
        
        elif type == 'vcf':
            with open(file) as f:
                print('vcf parsing is not yet implemented, coming soon...')
                sys.exit()
                # implement vcf object? or generate list

    except IOError as e:
        errno, strerror = e.args
        logging.fatal(e)
        print('Unable to open file, Error {}: {}'.format(errno,strerror))

# # # # takes python object from clinvar and puts info into dataframe
def parse_to_tbl(query_results):
    # build dataframe
    for result in query_results:
        # grab easy data
        parse_submissions(submission_subtree)
    return

# # # # takes subtree from main for submitters and calculates submission stats
def parse_submissions(submission_subtree):
    # build dataframe
    # aggregate functions for data
    return

# # # # outfile generation with vcf option
def output_files(type, ):
    # tbl output
    if type == 'vcf':
        vcf.output_vcf()
    return
    
# # # # 
def main():
    logging.basicConfig(level=logging.DEBUG)
    
    args = getargs()
    Entrez.email = args.email
    logging.debug('CLI inputs are {} {} {}'.format(args.type, args.email, args.input))
    
    for file in args.input:
        vid_list = []
        query_results = []
        input_selection(args.type, file, vid_list)
        logging.debug('vid_list to query clinvar: total {} items'.format(len(vid_list)))
        getncbi.get_entrez_xml(vid_list, query_results)
        logging.debug('the keys for query_results[0]: {}'.format(query_results[0].keys()))
    

if __name__ == '__main__':
    main()