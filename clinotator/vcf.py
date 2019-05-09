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
import pandas as pd
import global_vars as g


# error logging function plus config
def error_handling():
    return ' {}. {}, line: {}'.format(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2].tb_lineno)

# generates header for outfile, also counts input info and header lines
def parse_header(file_object, outprefix):

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

# process input vcf file, return rsids for query and vcf_tbl for output
def vcf_prep(file_object, outprefix):
    header_count = parse_header(file_object, outprefix)
    vcf_tbl = pd.read_csv(file_object, sep='\t', skiprows=header_count,
                          dtype=str)
    logging.debug('vcf_tbl shape -> {}'.format(vcf_tbl.shape))
    vcf_list = vcf_tbl.ID.values[vcf_tbl.ID.values != '.'].tolist()
    return vcf_list, vcf_tbl

# concatenates new info into variants that match rsid and alt
def cat_info_column(info, rsid, alt, out_tbl):
    rsid_match = rsid.lstrip('rs')
    alt_list = alt.split(",")
    info_columns = ['VID', 'CVVT', 'CVAL', 'CVCS', 'CVSZ', 'CVNA', 'CVDS',
                    'CVLE', 'CTRS', 'CTAA', 'CTPS', 'CTRR']
    logging.debug('rsid: {} alt_list: {}'.format(rsid_match, alt_list))
    # logging.debug('out_tbl shape -> {}'.format(out_tbl.shape))
    info_tbl = out_tbl.loc[(out_tbl['rsID'].astype('str') == rsid_match)
                           & out_tbl['CVAL'].isin(alt_list)].copy()
    
    if len(info_tbl.index) > 0:
        info_tbl.replace({'CVCS': {',': '%2C', ';': '%3B'},
                          'CVDS': {',': '%2C', ';': '%3B'}},
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
