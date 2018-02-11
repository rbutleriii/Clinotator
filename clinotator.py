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

import argparse, re, logging, sys
import pandas as pd
from Bio import Entrez

__version__ = "0.1.0"

# error logging function plus config
def error_handling():
    # import logging
    # import sys
    # logging.basicConfig(level=logging.WARNING)
    
    return ' {}. {}, line: {}'.format(sys.exc_info()[0],
                                           sys.exc_info()[1],
                                           sys.exc_info()[2].tb_lineno)
# argparse function
def getargs():
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

def input_selection(type, file, vid_list):
    try:
        if type == 'vid':
            with open(file) as f:
                vid_list = [line.rstrip('\n') for line in f]
            return vid_list
        
        elif type == 'rsid':
            with open(file) as f:
                rsid_list = [line.rstrip('\n') for line in f]
                rsid_to_vid(rsid_list, vid_list)
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

# looking up rsids and converting them to vids
def rsid_to_vid(rsid_list, vid_list):
    # import Bio.Entrez as Entrez
    # Entrez.email = args.email
    # Entrez.tool = 'Clinotator' # preferred by NCBI

    logging.debug('rsid_list {}'.format(rsid_list))
    query = Entrez.read(Entrez.epost('snp', id=",".join(rsid_list)))
    webenv = query['WebEnv']
    query_key = query['QueryKey']
    logging.debug('returned webenv: {} returned query key {}'.format(webenv, query_key))
    query_elink = Entrez.read(Entrez.elink(dbfrom='snp', db='clinvar', webenv=webenv, query_key=query_key)) # REPLACE WITH WEB HISTORY
    for link_list in query_elink:
        for linksetdb in link_list['LinkSetDb']:
            for link_id in linksetdb['Link']:
                vid_list.append(link_id['Id'])
                # each link query returns an 'IdList' of one, so always link_list['IdList'][0].
                # May not work with webenv, as join instead returns one search result for all
                logging.debug('rsid {} returns VID: {}'.format(link_list['IdList'][0], link_id['Id']))
    logging.debug('vid_list -> {}'.format(vid_list))
    
    return vid_list
    
def main():
    Entrez.tool = 'Clinotator' # preferred by NCBI
    logging.basicConfig(level=logging.DEBUG)
    
    args = getargs()
    Entrez.email = args.email
    logging.debug('CLI inputs are {} {} {}'.format(args.type, args.email, args.input))
    
    for file in args.input:
        vid_list = []
        input_selection(args.type, file, vid_list)
        logging.debug('vid_list to query clinvar: total {} items'.format(len(vid_list)))
    
    

if __name__ == '__main__':
    main()