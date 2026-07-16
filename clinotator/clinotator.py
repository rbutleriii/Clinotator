#!/usr/bin/env python

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
from typing import List, Union
import pandas as pd
import numpy as np
import Bio.Entrez as Entrez
import getncbi
import vcf
import variation
import global_vars as g


__version__ = g.__version__

def error_handling() -> str:
    """Format the currently-handled exception for logging.

    Must be called from within an ``except`` block, since it reads the
    exception info of the exception currently being handled.

    Returns:
        A message with the exception type, message, and line number.
    """
    return ' {}. {}, line: {}'.format(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2].tb_lineno)

def getargs() -> argparse.Namespace:
    """Define and parse the command-line arguments for clinotator.py.

    Returns:
        The parsed arguments, including ``input``, ``type``, ``email``,
        ``outprefix``, ``log``, and ``long_log``.
    """
    parser = argparse \
            .ArgumentParser(prog='clinotator.py',
                            formatter_class=argparse.RawTextHelpFormatter,
                            description='Clinical interpretation of ambiguous'
                                        ' ClinVar annotations')
    parser.add_argument('--log', action='store_true', help='create logfile')
    parser.add_argument('--long-log', action='store_true',
                        help='create detailed logfile')
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

def log_opts(log: bool, long_log: bool, outprefix: str) -> None:
    """Configure the logging module based on the chosen CLI log options.

    Args:
        log: If True, write INFO-level logs to ``{outprefix}.log``.
        long_log: If True, write DEBUG-level logs to ``{outprefix}.log``.
            Takes effect only when ``log`` is False.
        outprefix: Prefix used to build the log file name.
    """
    if log:
        logging.basicConfig(level=logging.INFO,
                            filename="{}.log".format(outprefix))
    elif long_log:
        logging.basicConfig(level=logging.DEBUG,
                            filename="{}.log".format(outprefix))
    else:
        logging.basicConfig(level=logging.INFO)

def input_selection(
        file_type: str,
        file: str,
        outprefix: str,
        query_results: List[str]) -> Union[pd.DataFrame, bool]:
    """Read an input file and dispatch it to the right NCBI query routine.

    Args:
        file_type: One of ``'vid'``, ``'rsid'``, or ``'vcf'``.
        file: Path to the input file.
        outprefix: Prefix used for any files written during vcf processing.
        query_results: List to be populated in-place with raw ClinVar XML
            responses fetched from NCBI.

    Returns:
        The parsed vcf table (``pandas.DataFrame``) if ``file_type`` is
        ``'vcf'``, otherwise ``False``. Returns ``None`` if the file could
        not be opened or another exception occurred.
    """
    try:
        with open(file) as f:

            if file_type == 'vid':
                id_list = np.unique([line.rstrip('\r\n') for line in f])
                getncbi.get_ncbi_xml(file_type, id_list, query_results)
                return False
        
            elif file_type == 'rsid':
                id_list = np.unique([line.lstrip('rsRS')
                                     .rstrip('\r\n') for line in f])
                getncbi.get_ncbi_xml(file_type, id_list, query_results)
                return False
            
            elif file_type == 'vcf':
                vcf_list, vcf_tbl = vcf.vcf_prep(f, outprefix)
                id_list = np.unique([item.lstrip('rs') for item in vcf_list])
                getncbi.get_ncbi_xml('rsid', id_list, query_results)
                return vcf_tbl
                
    except IOError as e:
        logging.critical(error_handling())
        print('Unable to open file, {}'.format(error_handling()))
    
    except:
        logging.critical(error_handling())

def explode(
        df: pd.DataFrame,
        lst_cols: Union[str, List[str]],
        fill_value: str = '') -> pd.DataFrame:
    """Explode list-valued cells in a DataFrame into separate rows.

    Function courtesy of @MaxU on stackoverflow question-12680754.

    Args:
        df: Source DataFrame containing one or more columns of lists.
        lst_cols: Column name, or list of column names, whose list values
            should be exploded in lockstep (all given columns must have
            lists of the same length per row).
        fill_value: Value used to fill columns for rows where the list in
            ``lst_cols[0]`` was empty.

    Returns:
        A new DataFrame with one row per list element, in the original
        column order.
    """
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

def output_files(
        vcf_tbl: Union[pd.DataFrame, bool],
        variant_objects: List['variation.VariationClass'],
        outprefix: str) -> None:
    """Write the per-variant results table, and an annotated vcf if given.

    Args:
        vcf_tbl: The vcf DataFrame produced by ``vcf.vcf_prep``, or ``False``
            if the input was not a vcf file.
        variant_objects: Parsed ``VariationClass`` instances, one per
            variant, used to build the output table.
        outprefix: Prefix used for the ``{outprefix}.tsv`` and, if
            applicable, ``{outprefix}.anno.vcf`` output files.
    """
    columnz = ['VID', 'CVVT', 'rsID', 'CVAL', 'vcfmatch', 'CVCS', 'CVSZ',
               'CVNA', 'CVDS', 'CVLE', 'CTRS', 'CTAA', 'CTPS', 'CTRR']

    result_tbl = pd.DataFrame([{fn: getattr(variant, fn) for fn in columnz}
            for variant in variant_objects])
    result_tbl = result_tbl[columnz]
    result_tbl['VID'] = result_tbl['VID'].astype(int)
    result_tbl.sort_values(by='VID', inplace=True)
    logging.debug('result_tbl shape -> {}'.format(result_tbl.shape))

    out_tbl = explode(result_tbl, ['rsID', 'CVAL'], fill_value='.')
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
    

def main() -> None:
    """Run the Clinotator CLI: parse args, query NCBI, and write results.

    Iterates over every input file given on the command line, queries
    NCBI/ClinVar for each, computes Clinotator metrics, and writes a
    ``.tsv`` (and, for vcf input, an annotated ``.anno.vcf``) per file.
    """
    args = getargs()
    log_opts(args.log, args.long_log, args.outprefix)
    logging.info('Run date: {}'
                 .format(datetime.date.today().strftime('%Y-%m-%d')))
    logging.debug('CLI inputs are {} {} {} {}'
                  .format(args.type, args.email, args.input, args.outprefix))
    Entrez.email = args.email
    
    for file in args.input:
        logging.info('Starting on {}'.format(file))
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
        logging.info('file written to {}.tsv'.format(outprefix))
    sys.exit()

if __name__ == '__main__':
    main()
