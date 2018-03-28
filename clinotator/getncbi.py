#!/usr/bin/env python

'''
Below are the functions for querying the entrez databases using biopython
Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import logging
import sys
import time
import Bio.Entrez as Entrez
import global_vars as g
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


Entrez.tool = g.etool

# timing compatibility with python2
def timing_tool():
    try:
        return time.perf_counter()
    except:
        return time.time()

# epost list to e-utilities history server for target db.
# Can't elink over 1000, no max on epost
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

    logging.debug('returned webenv: {} and query key: {}'
                  .format(webenv, query_key))
    return webenv, query_key
    
# Utilized for batch efetch/esummary from ncbi using history server.
def batch_ncbi(query_type, query_results, id_list, **kwargs):
    count = len(id_list)

    for start in range(0, count, g.efetch_batch):
        end = min(count, start + g.efetch_batch)
        logging.debug('{} run with {}'
                      .format(query_type, dict(retstart=start, **kwargs)))
        print("Going to download record {} to {}".format(start + 1, end))

        for attempt in range(4): # HTTPError and retry three times
            try:
                fetch_handle = getattr(Entrez, query_type) \
                        (**dict(retstart=start, **kwargs))
            except ValueError as oops:
                logging.warning('Empty set: Likely total number = batch size')
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Attempt {} of 4".format(attempt + 1))
                    logging.warning(err)
                    time.sleep(15)
                elif err.code == 400:
                    print("Attempt {} of 4".format(attempt + 1))
                    logging.warning(err)
                    sys.tracebacklimit = None
                    time.sleep(15)
                else:
                    raise
            else:
                break
        else:
            logging.critical('Failed to connect to NCBI 4 times exiting...')
            sys.exit()

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

# E-utilities batch queries w/o history server, useful for elink.
# Unreliable over batch size of 1000
def batch_nohist(query_type, id_list, **kwargs):
    count = len(id_list)
    result_list = []

    for start in range(0, count, g.elink_batch):
        end = min(count, start + g.elink_batch)
        logging.debug('{} run with {}'
                      .format(query_type,
                              dict(id=",".join(id_list[start:end]), **kwargs)))
        print("Looking up VIDs for rsIDs {} to {}".format(start + 1, end))
        
        for attempt in range(4): # HTTPError and retry three times
            try:
                fetch_handle = getattr(Entrez, query_type) \
                        (**dict(id=",".join(id_list[start:end]), **kwargs))
            except ValueError as oops:
                logging.warning('Empty set: Likely total number = batch size')
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Attempt {} of 4".format(attempt + 1))
                    logging.warning(err)
                    time.sleep(15)
                elif err.code == 400:
                    print("Attempt {} of 4".format(attempt + 1))
                    logging.warning(err)
                    sys.tracebacklimit = None
                    time.sleep(15)
                else:
                    raise
            else:
                break
        else:
            logging.critical('Failed to connect to NCBI 4 times, exiting...')
            sys.exit()

        record = Entrez.read(fetch_handle)
        result_list.extend(
                [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']])
    return result_list

# getting xml variation files for query_results list, 
# main
def get_ncbi_xml(file_type, id_list, query_results):
    logging.debug('{} list -> {}'.format(file_type, id_list))
    start = timing_tool()
    staging_list = []

    if file_type == 'rsid':
        id_list = list(filter(None, id_list)) # costs ~ 6 sec per 50k list
        vid_list = batch_nohist('elink', id_list, db='clinvar',
                                dbfrom='snp', linkname='snp_clinvar')
    elif file_type == 'vid':
        vid_list = list(filter(None, id_list)) # costs ~ 6 sec per 50k list
    else:
        logging.critical('Error: Incorrect file_type argument in get_rsid_xml ->'
                      ' {}'.format(file_type))

    webenv1, query_key1 = post_ncbi(file_type, 'epost', db='clinvar',
                                    id=",".join(vid_list))
    batch_ncbi('efetch', staging_list, vid_list, db='clinvar',
               rettype='variation', retmax=g.efetch_batch, webenv=webenv1,
               query_key=query_key1)
    [query_results.append(batch) for batch in staging_list] 
    stop = timing_tool()
    logging.info('Download time: {} min, Batches run -> {}'
                  .format(((stop-start)/60), len(query_results)))
    return
