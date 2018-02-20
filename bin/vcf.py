#!/usr/bin/env python3

'''
Below are the functions for manipulation of vcf files

Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import re
import pandas as pd

# generates header for outfile, also counts input header lines
def parse_header(file_object, outprefix):
    try:
        with open('{}.anno.vcf'.format(outprefix), 'w') as f:
            header = []
            info_list = []
            for index, line in enumerate(next(f) for x in range(200)):
                m = re.match('##(\w+)=', line)
                if m and m.group(1) == 'INFO':
                    header.append(line)
                    info_list.append(index)
                elif m and line.group(1) != 'INFO':
                    header.append(line)
                else:
                    continue
            header_count = len(header)
            header[max(info_list):max(info_list)] = new_headers
            f.write(header)
        return header_count
        
    except IOError as e:
        logging.fatal(error_handling())
        print('Unable to open file, {}'.format(error_handling()))
    
    except:
        logging.fatal(error_handling())

def body_to_tbl():
    pass

# setup for ncbi query, vcf outfile
def file_prep(file_object, outprefix):
    header_count = parse_header(file_object, outprefix)
    body_to_tbl(file_object, header_count)
    return vcf_list, vcf_tbl

# # # # writing the new data to the INFO track and out.vcf
def output_vcf():
    return

if __name__ == '__main__':
    with open('test.vcf', 'r') as file_object:
        parse_header(file_object, 'test_header')