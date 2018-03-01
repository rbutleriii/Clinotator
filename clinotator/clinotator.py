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


import os
import argparse
import logging
import sys
import datetime
import pandas as pd
import numpy as np
import Bio.Entrez as Entrez
import getncbi
import vcf
import variation
import global_vars as g


__version__ = g.__version__

# error logging function plus config
def error_handling():
    return ' {}. {}, line: {}'.format(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2].tb_lineno)

# argparse function
def getargs():
    parser = argparse \
            .ArgumentParser(prog='clinotator.py',
                            formatter_class=argparse.RawTextHelpFormatter,
                            description='Clinical interpretation of ambiguous'
                                        ' ClinVar annotations')
    parser.add_argument('--log', action='store_true', help='create logfile')
    parser.add_argument('-o', metavar='prefix', dest='outprefix',
                        default='clinotator',
                        help='choose an alternate prefix for outfiles')
    parser.add_argument("--version", action='version',
                        version='\n'.join(['Clinotator v'
                                          + __version__, __doc__]))
    parser.add_argument("input", metavar=('file'), nargs='+',
                        help="input file(s) (returns outfile for each)")
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed \
            .add_argument('-e', dest='email', required=True,
                          help='NCBI requires an email for database queries')
    requiredNamed.add_argument('-t', dest='type', required=True,
                               choices=['vid', 'rsid', 'vcf'],
                               help='vid - ClinVar Variation ID list\n'
                                    'rsid - dbSNP rsID list\n'
                                    'vcf - vcf file (output vcf generated)')
    return parser.parse_args()

# how to handle file types, returns vcf_tbl or False for output 
def input_selection(file_type, file, outprefix, query_results):
    try:
        with open(file) as f:

            if file_type == 'vid':
                id_list = [line.rstrip('\n') for line in f]
                getncbi.get_ncbi_xml(file_type, id_list, query_results)
                return False
        
            elif file_type == 'rsid':
                id_list = [line.lstrip('rsRS').rstrip('\n') for line in f]
                getncbi.get_ncbi_xml(file_type, id_list, query_results)
                return False
            
            elif file_type == 'vcf':
                vcf_list, vcf_tbl = vcf.vcf_prep(f, outprefix)
                id_list = [item.lstrip('rs') for item in vcf_list]
                getncbi.get_ncbi_xml('rsid', id_list, query_results)
                return vcf_tbl
                
    except IOError as e:
        logging.fatal(error_handling())
        print('Unable to open file, {}'.format(error_handling()))
    
    except:
        logging.fatal(error_handling())

# exploding list cells in dataframe into separate rows.
def explode(df, lst_cols, fill_value=''):
    # make sure `lst_cols` is a list
    if lst_cols and not isinstance(lst_cols, list):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)

    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()

    if (lens > 0).all():
        # ALL lists in cells aren't empty
        return pd.DataFrame({
            col:np.repeat(df[col].values, df[lst_cols[0]].str.len())
            for col in idx_cols
        }).assign(**{col:np.concatenate(df[col].values) for col in lst_cols}) \
          .loc[:, df.columns]
    else:
        # at least one list in cells is empty
        return pd.DataFrame({
            col:np.repeat(df[col].values, df[lst_cols[0]].str.len())
            for col in idx_cols
        }).assign(**{col:np.concatenate(df[col].values) for col in lst_cols}) \
          .append(df.loc[lens==0, idx_cols]).fillna(fill_value) \
          .loc[:, df.columns]

# outfile generation with vcf option
def output_files(vcf_tbl, variant_objects, outprefix):
    columnz = ['VID', 'CVVT', 'RSID', 'CVMA', 'vcfmatch', 'CVCS', 'CVSZ',
               'CVNA', 'CVDS', 'CVLE', 'CTRS', 'CTAA', 'CTPS', 'CTRR']

    result_tbl = pd.DataFrame([{fn: getattr(variant, fn) for fn in columnz}
        for variant in variant_objects])
    result_tbl = result_tbl[columnz]
    # result_tbl = result_tbl.sort_values(by='VID', axis=1)
    logging.debug('result_tbl shape -> {}'.format(result_tbl.shape))

    out_tbl = explode(result_tbl, ['RSID', 'CVMA'], fill_value='.')
    out_tbl.to_csv('{}.tsv'.format(outprefix), sep='\t', na_rep='.',
                   index=False)
    logging.debug('out_tbl shape -> {}'.format(out_tbl.shape))

    if isinstance(vcf_tbl, pd.DataFrame):
        vcf_tbl['INFO'] = vcf_tbl \
            .apply(lambda x: vcf.cat_info_column(x['INFO'], x['ID'], x['ALT'],
                                                 out_tbl), axis=1)
        vcf_tbl.to_csv('{}.anno.vcf'.format(outprefix), sep='\t', mode='a',
                       index=False, na_rep='.')
    return
    
def main():
    args = getargs()
    
    if args.log:
        logging.basicConfig(level=logging.DEBUG, filename="clinotator.log")
    else:
        logging.basicConfig(level=logging.WARN)

    Entrez.email = args.email
    logging.debug('CLI inputs are {} {} {} {}'
                  .format(args.type, args.email, args.input, args.outprefix))
    
    for file in args.input:
        query_results = []
        variant_objects = []
        base = os.path.basename(file)
        outprefix = args.outprefix + '.' + os.path.splitext(base)[0]

        vcf_tbl = input_selection(args.type, file, outprefix, query_results)
        logging.debug('the total # of query_results: {}'
                      .format(len(query_results)))
        
        variation.query_parsing(variant_objects, query_results)
        logging.debug('the total # of variant_objects: {}'
                      .format(len(variant_objects)))
        
        output_files(vcf_tbl, variant_objects, outprefix)
        logging.debug('file written to {}.tsv'.format(outprefix))
        sys.exit()

if __name__ == '__main__':
    main()