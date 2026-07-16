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
from typing import Any, List, Tuple
import Bio.Entrez as Entrez
import global_vars as g
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


Entrez.tool = g.etool

def timing_tool() -> float:
    """Get a monotonic timestamp, compatible with Python 2 and 3.

    Returns:
        A timestamp in fractional seconds from ``time.perf_counter`` (or
        ``time.time`` as a Python 2 fallback). Only meaningful relative to
        another call to this function.
    """
    try:
        return time.perf_counter()
    except:
        return time.time()

def http_attempts(query_type: str, **kwargs: Any) -> Any:
    """Call a Bio.Entrez function, retrying on transient HTTP errors.

    Args:
        query_type: Name of the ``Bio.Entrez`` function to call (e.g.
            ``'epost'``, ``'efetch'``, ``'elink'``).
        **kwargs: Keyword arguments forwarded to the ``Bio.Entrez``
            function.

    Returns:
        The open handle returned by the successful ``Bio.Entrez`` call.

    Raises:
        HTTPError: Re-raised for HTTP error codes other than 400, 429, or
            5xx, which are retried instead.

    Note:
        Exits the process via ``sys.exit()`` if all 4 attempts fail.
    """
    for attempt in range(4): # HTTPError and retry three times
        try:
            fetch_handle = getattr(Entrez, query_type)(**kwargs)
        except ValueError as oops:
            logging.warning('Empty set: Likely total number = batch size')
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                logging.info("Attempt {} of 4".format(attempt + 1))
                logging.warning(err)
                time.sleep(15)
            elif err.code == 400:
                logging.info("Attempt {} of 4".format(attempt + 1))
                logging.warning(err)
                sys.tracebacklimit = None
                time.sleep(15)
            elif err.code == 429:
                logging.info("Attempt {} of 4".format(attempt + 1))
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
    return fetch_handle

def post_ncbi(
        file_type: str,
        query_type: str,
        **kwargs: Any) -> Tuple[str, str]:
    """Submit a query to NCBI and retrieve its history-server reference.

    Args:
        file_type: One of ``'vid'``, ``'rsid'``, or ``'vcf'``; used only to
            select how the response is unpacked.
        query_type: Name of the ``Bio.Entrez`` function to call, e.g.
            ``'epost'`` or ``'elink'``.
        **kwargs: Keyword arguments forwarded to the ``Bio.Entrez``
            function.

    Returns:
        A ``(webenv, query_key)`` tuple identifying the posted query on
        NCBI's history server, for use in later batched efetch/elink calls.
    """
    logging.debug('{} to ncbi using kwargs: {}'.format(query_type, kwargs))
    fetch_handle = http_attempts(query_type, **kwargs)
    query = Entrez.read(fetch_handle)
    fetch_handle.close()

    if file_type == 'rsid' and query_type == 'elink':
        webenv = query[0]['WebEnv']
        query_key = query[0]['LinkSetDbHistory'][0]['QueryKey']
    else:    
        webenv = query['WebEnv']
        query_key = query['QueryKey']

    logging.debug('returned webenv: {} and query key: {}'
                  .format(webenv, query_key))
    return webenv, query_key
    
def batch_ncbi(
        query_type: str,
        query_results: List[str],
        id_list: List[str],
        **kwargs: Any) -> None:
    """Fetch records from NCBI in batches, using the history server.

    Args:
        query_type: Name of the ``Bio.Entrez`` function to call, typically
            ``'efetch'``.
        query_results: List to be populated in-place with each batch's raw
            response text.
        id_list: The full list of IDs being fetched; only its length is
            used, to compute batch boundaries.
        **kwargs: Keyword arguments forwarded to the ``Bio.Entrez`` function
            (e.g. ``db``, ``rettype``, ``webenv``, ``query_key``).
    """
    count = len(id_list)

    time.sleep(0.35)
    for start in range(0, count, g.efetch_batch):
        end = min(count, start + g.efetch_batch)
        logging.info("Going to download record {} to {}"
                     .format(start + 1, end))

        logging.debug('{} run with {}'
                      .format(query_type, dict(retstart=start, **kwargs)))
        fetch_handle = http_attempts(query_type,
                                     **dict(retstart=start, **kwargs))
        # Ideally Entrez.read would import as python object, and validate with
        # ClinVar, but ClinVar currently deoesn't declare DOCTYPE for 
        # validation in their xml header. They have .xsd for this purpose, but
        # Entrez.read can't use this w/o a declaration to validate. Something
        # of a security flaw, but hopefully NCBI is in control of their xml
        # and it is injection free (hopefully...). Ultimately, the best option
        # for now is to utilize ElementTree for parsing. >Revisit<
        # data = Entrez.read(fetch_handle)
        # fetch_handle.close()
        query_results.append(fetch_handle.read())
        fetch_handle.close()
    return

def batch_local(
        file_type: str,
        query_type: str,
        id_list: List[str],
        **kwargs: Any) -> List[str]:
    """Locally batch an ID list, then epost each batch before querying.

    ``elink`` can't handle a batch size over 1000 and doesn't support
    ``retstart``/``retmax``, so this is the only way to run ``elink``,
    ``egquery``, etc. over a large ID list: batch locally, then
    epost -> query for each batch.

    Args:
        file_type: One of ``'vid'``, ``'rsid'``, or ``'vcf'``; forwarded to
            ``post_ncbi`` to control response unpacking.
        query_type: Name of the ``Bio.Entrez`` function to call after
            posting each batch, typically ``'elink'``.
        id_list: Full list of IDs (e.g. rsIDs) to look up.
        **kwargs: Keyword arguments forwarded to ``post_ncbi`` (for the
            epost call) and to ``http_attempts`` (for the query call),
            e.g. ``dbfrom``, ``db``, ``linkname``.

    Returns:
        The combined list of linked IDs (e.g. ClinVar VIDs) found for all
        batches of ``id_list``.
    """
    count = len(id_list)
    logging.debug('total rsIDs to post: {}'.format(count))
    result_list = []

    for start in range(0, count, g.elink_batch):
        end = min(count, start + g.elink_batch)
        logging.info("Looking up VIDs for rsIDs {} to {}"
                     .format(start + 1, end))
        webenv1, query_key1 = post_ncbi(
                file_type, 'epost', db=kwargs['dbfrom'],
                id=",".join(id_list[start:end]))
        logging.debug('{} run with {}'.format(query_type, dict(
                webenv=webenv1, query_key=query_key1, **kwargs)))
        fetch_handle = http_attempts(query_type, **dict(
                webenv=webenv1, query_key=query_key1, **kwargs))
        record = Entrez.read(fetch_handle)
        fetch_handle.close()
        try:
            result_list.extend(
                    [lnk['Id'] for lnk in record[0]['LinkSetDb'][0]['Link']])
        except IndexError:
            logging.info('No VIDs for rsIDs {} to {}'.format(start + 1, end))
            pass
        continue
        logging.debug('length result list: {}'.format(len(result_list)))
        time.sleep(0.35)
    return result_list

def get_ncbi_xml(
        file_type: str,
        id_list: List[str],
        query_results: List[str]) -> None:
    """Resolve IDs to ClinVar VIDs (if needed) and fetch their XML records.

    This is the main entry point of this module: given a list of vid or
    rsid identifiers, it resolves them to ClinVar Variation IDs if
    necessary, then downloads the corresponding ClinVar variation-report
    XML in batches.

    Args:
        file_type: One of ``'rsid'`` or ``'vid'``.
        id_list: The list of identifiers to resolve/fetch.
        query_results: List to be populated in-place with raw ClinVar XML
            responses fetched from NCBI.
    """
    logging.debug('{} list -> {}'.format(file_type, id_list))
    start = timing_tool()
    staging_list = []

    if file_type == 'rsid':
        id_list = list(filter(None, id_list)) # costs ~ 6 sec per 50k list
        vid_list = batch_local(file_type, 'elink', id_list, db='clinvar',
                                dbfrom='snp', linkname='snp_clinvar')
    elif file_type == 'vid':
        vid_list = list(filter(None, id_list)) # costs ~ 6 sec per 50k list
    else:
        logging.critical('Error: Incorrect file_type argument in get_rsid_xml'
                ' -> {}'.format(file_type))

    webenv1, query_key1 = post_ncbi(file_type, 'epost', db='clinvar',
                                    id=",".join(vid_list))
    batch_ncbi('efetch', staging_list, vid_list, db='clinvar',
               rettype='vcv', is_variationid='true', retmax=g.efetch_batch,
               webenv=webenv1, query_key=query_key1)
    [query_results.append(batch) for batch in staging_list] 
    stop = timing_tool()
    logging.info('Download time: {} min, Batches run -> {}'
                 .format(((stop-start)/60), len(query_results)))
    return
