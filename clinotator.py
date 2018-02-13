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


import argparse
import logging
import sys
import pandas as pd
import bin.getncbi as getncbi
import bin.vcf as vcf
import bin.variation as variation
from Bio import Entrez


__version__ = "0.0.1"

# error logging function plus config
def error_handling():
    return ' {}. {}, line: {}'.format(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2].tb_lineno)

# argparse function
def getargs():
    parser = argparse.ArgumentParser(prog='clinotator',
                                     formatter_class=
                                     argparse.RawTextHelpFormatter,
                                     description='Clinical interpretation ' \
                                     'of ambiguous ClinVar annotations')
    parser.add_argument("--version", action='version',
                        version='\n'.join(['Clinotator v'
                                          + __version__, __doc__]))
    parser.add_argument('-o', dest='outprefix', default='clinotator.file',
                        nargs='?', const='clinotator.file',
                        help='Choose an alternate prefix for outfiles')
    parser.add_argument("input", metavar=('file'), nargs='+',
                        help="input file(s) (returns outfile for each)")
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-t', dest='type', required=True,
                               choices=['vid', 'rsid', 'vcf'],
                               help='vid - ClinVar Variation ID list\n' \
                                    'rsid - dbSNP rsID list\n' \
                                    'vcf - vcf file (output vcf generated)')
    requiredNamed.add_argument('-e', dest='email', required=True,
                               help='NCBI requires an email for querying '\
                                    'their databases')
    return parser.parse_args()

# # # # how to handle file types 
def input_selection(file_type, file, query_results):
    try:
        with open(file) as f:
            id_list = [line.lstrip('rsRS').rstrip('\n') for line in f]

            if file_type == 'vid':
                getncbi.get_ncbi_xml(file_type, id_list, query_results)
        
            elif file_type == 'rsid':
                getncbi.get_ncbi_xml(file_type, id_list, query_results)
        
            elif file_type == 'vcf':
                print('vcf parsing is not yet implemented, coming soon...')
                sys.exit()
                # implement vcf object? or generate list
        return

    except IOError as e:
        errno, strerror = e.args
        logging.fatal(e)
        print('Unable to open file, Error {}: {}'.format(errno,strerror))
    
    except:
        logging.fatal(error_handling())

# # # # outfile generation with vcf option
def output_files(file_type, variant_objects, outprefix):
    columns = ['VID', 'RSID', 'CTCS', 'CTSZ', 'CTNA', 'CTDS', 'CTLU',
               'CTS', 'CTAA', 'CTWS', 'CTRR']
    result_tbl = pd.DataFrame([{fn: getattr(variant, fn) for fn in columns}
        for variant in variant_objects])
    result_tbl.sort_values(['VID', 'RSID'], inplace=True)
    result_tbl.drop_duplicates(inplace=True)
    result_tbl.to_csv('{}.tsv'.format(outprefix), sep='\t', na_rep='.',
                      index=False, tupleize_cols=True)
    if file_type == 'vcf':
        vcf.output_vcf()
    return
    
def main():
    logging.basicConfig(level=logging.DEBUG)
    args = getargs()
    Entrez.email = args.email
    logging.debug('CLI inputs are {} {} {} {}'
                  .format(args.type, args.email, args.input, args.outprefix))
    
    for file in args.input:
        query_results = []
        variant_objects = []
        input_selection(args.type, file, query_results)
        logging.debug('the total # of query_results: {}'
                      .format(len(query_results)))
        
        for batch in query_results:
            clinvarresult = ET.fromstring(batch)
            for variationreport in clinvarresult.iter('VariationReport'):
                variant = VariationClass(variationreport)
                variant_objects.append(variant)
                # logging.debug('something')
        
        outprefix = args.outprefix
        output_files(args.type, variant_objects, outprefix)
        # logging.debug('something')
        sys.exit()

if __name__ == '__main__':
    main()