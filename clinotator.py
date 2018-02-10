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
import Bio.Entrez as Entz

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


# looking up rsids and converting them to vids
def rsid_to_vid(rsid_list):
    # import Bio.Entrez as Entz
    # Entz.email = args.email
    # Entz.tool = 'Clinotator' # preferred by NCBI

    results = Entz.read(Entz.epost('snp' , id = ','.join(rsid_list)))
    webenv = results['WebEnv']
    query_key = results['QueryKey']
    
    return
    
def main():
    logging.basicConfig(level=logging.DEBUG)
    args = getargs()
    logging.debug('CLI inputs are {} {} {}'.format(args.type, args.email, args.input))
    Entz.tool = 'Clinotator' # preferred by NCBI
    Entz.email = args.email

if __name__ == '__main__':
    main()