#!/usr/bin/env python

'''
Below are the functions for manipulation of vcf files

Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import re
import logging
import sys
import datetime
from typing import IO, List, Tuple
import pandas as pd
import global_vars as g


def error_handling() -> str:
    """Format the currently-handled exception for logging.

    Must be called from within an ``except`` block, since it reads the
    exception info of the exception currently being handled.

    Returns:
        A message with the exception type, message, and line number.
    """
    return ' {}. {}, line: {}'.format(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2].tb_lineno)

def parse_header(file_object: IO[str], outprefix: str) -> int:
    """Copy and extend a vcf header into the annotated output vcf.

    Reads header lines from the input vcf, inserts Clinotator's INFO field
    definitions just after the existing INFO definitions, adds a
    Clinotator run-metadata line, and writes the result to
    ``{outprefix}.anno.vcf``.

    Args:
        file_object: Open file handle for the input vcf, positioned at the
            start of the file. Its read cursor is reset to 0 before this
            function returns.
        outprefix: Prefix used to build the ``{outprefix}.anno.vcf`` output
            file name.

    Returns:
        The number of header lines read from the input file (used by the
        caller to know how many rows to skip when reading the body with
        pandas). Returns ``None`` if the file could not be opened or
        another exception occurred.
    """
    try:
        with open('{}.anno.vcf'.format(outprefix), 'w') as outfile:
            header = []
            info_list = []
            for index, line in zip(range(g.max_vcf_header_size), file_object):
                m = re.match('##([\w\-\.]+)=', line)
                
                if m and m.group(1) == 'INFO':
                    header.append(line)
                    info_list.append(index)
                elif m and m.group(1) != 'INFO':
                    header.append(line)
                else:
                    continue
            
            header_count = len(header)
            header[max(info_list) + 1 : max(info_list) + 1] = g.new_headers
            rundate = datetime.date.today().strftime('%Y-%m-%d')
            meta_clin_line = ('##annotation=CLINOTATORv{}_run_{}\n'
                              .format(g.__version__, rundate))
            header.insert(min(info_list), meta_clin_line)
            outfile.writelines(('{}'.format(item) for item in header))
            file_object.seek(0) # reset read cursor to beginning
            logging.debug('Initial header: {} INFO lines: {} Final header: {}'
                          .format(header_count, len(info_list), len(header)))
        return header_count
        
    except IOError as e:
        logging.fatal(error_handling())
        print('Unable to open file, {}'.format(error_handling()))
    
    except:
        logging.fatal(error_handling())

def vcf_prep(
        file_object: IO[str],
        outprefix: str) -> Tuple[List[str], pd.DataFrame]:
    """Write the annotated-vcf header and load the vcf body into a table.

    Args:
        file_object: Open file handle for the input vcf, positioned at the
            start of the file.
        outprefix: Prefix used to build the ``{outprefix}.anno.vcf`` output
            file name (via ``parse_header``).

    Returns:
        A ``(vcf_list, vcf_tbl)`` tuple: ``vcf_list`` is the list of
        non-missing variant IDs (rsIDs) from the vcf ``ID`` column, and
        ``vcf_tbl`` is the full vcf body as a DataFrame of strings.
    """
    header_count = parse_header(file_object, outprefix)
    vcf_tbl = pd.read_csv(file_object, sep='\t', skiprows=header_count,
                          dtype=str)
    logging.debug('vcf_tbl shape -> {}'.format(vcf_tbl.shape))
    vcf_list = vcf_tbl.ID.values[vcf_tbl.ID.values != '.'].tolist()
    return vcf_list, vcf_tbl

def cat_info_column(
        info: str,
        rsid: str,
        alt: str,
        out_tbl: pd.DataFrame) -> str:
    """Append Clinotator annotations to a vcf record's INFO field.

    Looks up the given rsID/alt-allele combination in the Clinotator
    results table and, if found, appends the corresponding Clinotator
    fields to the existing INFO string.

    Args:
        info: The existing INFO field string for this vcf record.
        rsid: The record's variant ID (e.g. ``'rs12345'``), as found in the
            vcf ``ID`` column.
        alt: The record's comma-separated ALT allele string.
        out_tbl: The Clinotator results table (as produced by
            ``clinotator.output_files``/``clinotator.explode``) to search
            for a matching rsID/allele.

    Returns:
        The INFO string with Clinotator fields appended if a match was
        found, otherwise the original ``info`` string unchanged.
    """
    rsid_match = rsid.lstrip('rs')
    alt_list = alt.split(",")
    info_tbl = pd.DataFrame()
    info_columns = ['VID', 'CVVT', 'CVAL', 'CVCS', 'CVSZ', 'CVNA', 'CVDS',
                    'CVLE', 'CTRS', 'CTAA', 'CTPS', 'CTRR']
    logging.debug('rsid: {} alt_list: {}'.format(rsid_match, alt_list))
    # logging.debug('out_tbl shape -> {}'.format(out_tbl.shape))
    if rsid_match != ".":
        info_tbl = out_tbl.loc[(out_tbl['rsID'].astype('str') == rsid_match)
                               & out_tbl['CVAL'].isin(alt_list)].copy()
    
    if len(info_tbl.index) > 0:
        info_tbl.replace({'CVCS': {',': '%2C', ';': '%3B', ' ': '_'},
                          'CVDS': {',': '%2C', ';': '%3B'},
                          'CVVT': {' ': '_'},
                          'CTPS': {' ': '_'}},
                         regex=True, inplace=True)
        new_info = ['{}={}'.format(x, info_tbl[x]
                    .to_csv(header=False, index=False, na_rep='.')
                    .strip('\n')) for x in info_columns]
        new_info = [string.replace('\n', ',') for string in new_info]
        logging.debug('{} had a match: {}'.format(rsid, new_info))
        info_list =';'.join([info] + new_info)
        logging.debug('adding {}'.format(info_list))
        return info_list
    else:
        logging.debug('keeping as is {}'.format(info))
        return info
    
# test
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logging.debug('this is a module test')
    with open('../test/test.vcf', 'r') as file_object:
        vcf_list, vcf_tbl = vcf_prep(file_object, 'test_header')
        logging.debug('vcf_list -> {}'.format(vcf_list))
        sample_tbl = pd.read_csv('../test/test.tbl', sep='\t', dtype=str)
        info_list = cat_info_column('NS=3;DP=11;AF=0.017', 'rs34376836', 'A',
                                    sample_tbl)
        info_list = cat_info_column('NS=3;DP=11;AF=0.017', '.', 'A',
                                    sample_tbl)
        info_list = cat_info_column('NS=3;DP=11;AF=0.017', 'rs118161496', 'C',
                                    sample_tbl)
        info_list = cat_info_column('NS=3;DP=11;AF=0.017', 'rs200401432', 'A',
                                    sample_tbl)
        info_list = cat_info_column('NS=3;DP=11;AF=0.017', 'rs16992990', 'A',
                                    sample_tbl)
