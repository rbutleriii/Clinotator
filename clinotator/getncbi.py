#!/usr/bin/env python3

'''
Below are the functions for querying the entrez databases using biopython
Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import logging
import sys
import Bio.Entrez as Entrez
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


Entrez.tool = 'Clinotator' # preferred by NCBI
batch_size = 4500

# getting xml variation files for query_results list, 
def get_ncbi_xml(file_type, id_list, query_results):
    logging.debug('{} list -> {}'.format(file_type, id_list))
    staging_list = []

    if file_type == 'rsid':
        webenv2, query_key2 = post_ncbi(file_type, 'epost', db='snp',
                                        id=",".join(id_list))
        webenv1, query_key1 = post_ncbi(file_type, 'elink', db='clinvar',
                                        webenv=webenv2, query_key=query_key2,
                                        dbfrom='snp', linkname='snp_clinvar',
                                        cmd='neighbor_history')
    elif file_type == 'vid':
        webenv1, query_key1 = post_ncbi(file_type, 'epost', db='clinvar',
                                        id=",".join(id_list))
    else:
        logging.fatal('Error: Incorrect file_type argument in get_rsid_xml ' \
                      '-> {}'.format(file_type))

    batch_ncbi('efetch', staging_list, id_list, db='clinvar',
               rettype='variation', retmax=batch_size, webenv=webenv1,
               query_key=query_key1)
    [query_results.append(batch) for batch in staging_list] 
    logging.debug('batches run -> {}'.format(len(query_results)))
    return

# epost or elink list to e-utilities history server for target db. 
def post_ncbi(file_type, query_type, **kwargs):
    logging.debug('{} to ncbi using kwargs: {}'.format(query_type, kwargs))
    handle = getattr(Entrez, query_type)(**kwargs)
    query = Entrez.read(handle)
    logging.debug(query)
    handle.close()

    if file_type == 'rsid' and query_type == 'elink':
        webenv = query[0]['WebEnv']
        query_key = query[0]['LinkSetDbHistory'][0]['QueryKey']
    else:    
        webenv = query['WebEnv']
        query_key = query['QueryKey']

    logging.debug('returned webenv: {} and query key: {}'.format(webenv,
                                                                 query_key))
    return webenv, query_key
    
# Utilized for batch efetching from ncbi. HTTPError and retry 
def batch_ncbi(query_type, query_results, id_list, **kwargs):
    count = len(id_list)

    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        logging.debug('{} run with {}'.format(query_type,
                                              dict(retstart=start, **kwargs)))
        print("Going to download record %i to %i" % (start+1, end))
        attempt = 0
        while attempt < 3:
            attempt += 1

            try:
                fetch_handle = getattr(Entrez,
                                       query_type)(**dict(retstart=start,
                                                          **kwargs))

            except ValueError as oops:
                logging.warning('Likely total = batch size')
                break

            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    logging.debug(err)
                    time.sleep(15)
                elif err.code == 400:
                    logging.warning(err)
                    sys.tracebacklimit = None
                    break
                else:
                    raise
        # Ideally Entrez.read would import as python object, and validate with
        # ClinVar, but ClinVar no longer declares DOCTYPE for validation in
        # their xml returns. They have .xsd for this purpose, but users
        # downloading their own validation file is silly, plus there is no
        # option with Entrez.read to specify the .xsd file. Ultimately, the
        # best option for now is to utilize ElementTree for parsing. >Revisit<
        # data = Entrez.read(fetch_handle)
        # print(data)
        # fetch_handle.close()
        query_results.append(fetch_handle.read())
    return
