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
                   '##INFO=<ID=CVVT,Number=1,Type=String,Description="ClinVar variant type">\n',
                   '##INFO=<ID=CVMA,Number=1,Type=String,Description="ClinVar minor allele">\n',
                   '##INFO=<ID=CVCS,Number=1,Type=String,Description="ClinVar clinical significance">\n',
                   '##INFO=<ID=CVSZ,Number=1,Type=Integer,Description="ClinVar stars">\n',
                   '##INFO=<ID=CVNA,Number=1,Type=Integer,Description="ClinVar number of clinical assertions">\n',
                   '##INFO=<ID=CVDS,Number=1,Type=String,Description="ClinVar conditions">\n',
                   '##INFO=<ID=CVLE,Number=1,Type=String,Description="ClinVar last evaluated">\n',
                   '##INFO=<ID=CTRS,Number=1,Type=Float,Description="Clinotator raw score">\n',
                   '##INFO=<ID=CTAA,Number=1,Type=Float,Description="Clinotator average clinical assertion age">\n',
                   '##INFO=<ID=CTWS,Number=1,Type=String,Description="Clinotator weighted significance">\n',
                   '##INFO=<ID=CTRR,Number=1,Type=String,Description="Clinotator reclassification recommendation">\n']
    try:
        with open('{}.anno.vcf'.format(outprefix), 'w') as outfile:
            header = []
            info_list = []
            for index, line in enumerate(next(file_object) for x in range(max_vcf_header_size)):
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
            logging.debug('Initial header: {} INFO lines: {} Final header: {}'.format(header_count, len(info_list), len(header)))
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
    vcf_list = vcf_tbl['ID'][~vcf_tbl['ID'].isin(['.'])].unique().tolist()
    return vcf_list, vcf_tbl

def cat_info_column(rsid, alt, info):
    pass
    
# # # # writing the new data to the INFO track and out.vcf
def output_vcf(out_tbl, vcf_tbl):
    return

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    with open('test.vcf', 'r') as file_object:
        vcf_list, vcf_tbl = vcf_prep(file_object, 'test_header')
        print(vcf_list)
        vcf_tbl