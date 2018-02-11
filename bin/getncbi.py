#!/usr/bin/env python3

'''
Below are the functions for querying the entrez databases using biopython
Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import logging
from Bio import Entrez
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2

Entrez.tool = 'Clinotator' # preferred by NCBI
batch_size = 3

# looking up rsids and converting them to vids
def rsid_to_vid(rsid_list, vid_list):
    logging.debug('rsid_list {}'.format(rsid_list))
    webenv1, query_key1 = post_ncbi(rsid_list, 'snp')
    staging_list = []
    batch_ncbi('elink', staging_list, rsid_list, 3, dbfrom='snp', db='clinvar', webenv=webenv1,
               query_key=query_key1, LinkName='snp_clinvar')
    [[vid_list.append(link['Id']) for link in batch['LinkSetDb'][0]['Link']]for batch in staging_list] 
    logging.debug('vid_list -> {}'.format(vid_list))
    
    return

# getting xml files for vid_list
def get_entrez_xml(vid_list, query_results):
    logging.debug('vid_list {}'.format(vid_list))
    webenv1, query_key1 = post_ncbi(vid_list, 'clinvar')
    staging_list = []
    batch_ncbi('efetch', staging_list, vid_list, 3, db='clinvar', rettype='variation', webenv=webenv1,
               query_key=query_key1)
    [query_results.append(batch) for batch in staging_list] 
    logging.debug('batches run -> {}'.format(len(query_results)))
    
    return

# Posts list to history for target db
def post_ncbi(post_list, db):
    query = Entrez.read(Entrez.epost(db, id=",".join(post_list)))
    webenv1 = query['WebEnv']
    query_key1 = query['QueryKey']
    logging.debug('returned webenv: {} returned query key: {}'.format(webenv1, query_key1))
    
    return webenv1, query_key1
    
# Utilized for batching from ncbi. Use history variables in kwargs. HTTPError and retry 
def batch_ncbi(query_type, return_list, post_list, batch_size, **kwargs):
    count = len(post_list)
    logging.debug('{} run with {}'.format(query_type, kwargs))
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                fetch_handle = getattr(Entrez, query_type)(**kwargs)
            except HTTPError as err:
                logging.debug(err)
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(15)
                else:
                    raise
        data = Entrez.read(fetch_handle)
        print(data)
        fetch_handle.close()
        return_list.append(data[0])
    
    return
