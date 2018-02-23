#!/usr/bin/env python3

'''
Below are the functions for manipulation of vcf files

Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import re
import logging
import sys
import pandas as pd


# error logging function plus config
def error_handling():
    return ' {}. {}, line: {}'.format(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2].tb_lineno)

# generates header for outfile, also counts input header lines
def parse_header(file_object, outprefix):
    max_vcf_header_size = 200
    new_headers = ['##INFO=<ID=VID,Number=1,Type=Integer,Description="ClinVar variation ID">\n',
                   '##INFO=<ID=CVVT,Number=A,Type=String,Description="ClinVar variant type">\n',
                   '##INFO=<ID=CVMA,Number=A,Type=String,Description="ClinVar minor allele">\n',
                   '##INFO=<ID=CVCS,Number=A,Type=String,Description="ClinVar clinical significance">\n',
                   '##INFO=<ID=CVSZ,Number=A,Type=Integer,Description="ClinVar stars">\n',
                   '##INFO=<ID=CVNA,Number=A,Type=Integer,Description="ClinVar number of clinical assertions">\n',
                   '##INFO=<ID=CVDS,Number=A,Type=String,Description="ClinVar conditions">\n',
                   '##INFO=<ID=CVLE,Number=A,Type=String,Description="ClinVar last evaluated">\n',
                   '##INFO=<ID=CTRS,Number=A,Type=Float,Description="Clinotator raw score">\n',
                   '##INFO=<ID=CTAA,Number=A,Type=Float,Description="Clinotator average clinical assertion age">\n',
                   '##INFO=<ID=CTWS,Number=A,Type=String,Description="Clinotator weighted significance">\n',
                   '##INFO=<ID=CTRR,Number=A,Type=String,Description="Clinotator reclassification recommendation">\n']
    try:
        with open('{}.anno.vcf'.format(outprefix), 'w') as outfile:
            header = []
            info_list = []
            for index, line in enumerate(next(file_object)
                                         for x in range(max_vcf_header_size)):
                m = re.match('##(\w+)=', line)
                
                if m and m.group(1) == 'INFO':
                    header.append(line)
                    info_list.append(index)
                elif m and m.group(1) != 'INFO':
                    header.append(line)
                else:
                    continue
            
            header_count = len(header)
            header[max(info_list) + 1 : max(info_list) + 1] = new_headers
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

# process input vcf file, return rsids for query and vcf_tbl for output
def vcf_prep(file_object, outprefix):
    header_count = parse_header(file_object, outprefix)
    vcf_tbl = pd.read_table(file_object, skiprows=header_count, dtype=str)
    logging.debug('vcf_tbl shape -> {}'.format(vcf_tbl.shape))
    vcf_list = vcf_tbl.loc[vcf_tbl['ID'].ne('.'), 'ID'].unique().tolist()
    return vcf_list, vcf_tbl

# concatenates new info into variants that match rsid and alt
def cat_info_column(info, rsid, alt, out_tbl):
    rsid_match = rsid.lstrip('rs')
    alt_list = alt.split(",")
    info_columns = ['VID', 'CVVT', 'CVMA', 'CVCS', 'CVSZ', 'CVNA', 'CVDS',
                    'CVLE', 'CTRS', 'CTAA', 'CTWS', 'CTRR']
    logging.debug('rsid: {} alt_list: {}'.format(rsid_match, alt_list))
    # logging.debug('out_tbl shape -> {}'.format(out_tbl.shape))
    info_tbl = out_tbl.loc[(out_tbl['RSID'].astype('str') == rsid_match)
                           & out_tbl['CVMA'].isin(alt_list)]
    
    if len(info_tbl.index) > 0:
        new_info = ['{}={}'.format(x, info_tbl[x]
                    .to_csv(header=None, index=False)
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
        sample_tbl = pd.read_table('../test/test.tbl', dtype=str)
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