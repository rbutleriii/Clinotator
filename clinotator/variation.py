#!/usr/bin/env python3

'''
Below are the functions for the Variation file object which contains
variant info for metrics.

Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import pprint
import logging
import datetime
import warnings
import numpy as np
import xml.etree.ElementTree as ET
import global_vars as g


# work through batch queries from ncbi and build list of variant objects
def query_parsing(variant_objects, query_results):
    
    for batch_index, batch in enumerate(query_results):
        clinvarresult = ET.fromstring(batch)
        for var_index, variationreport in enumerate(clinvarresult):
            variant = VariationClass(variationreport)
            variant_objects.append(variant)
        logging.debug('{} variant objects parsed in batch {}'
                      .format(var_index + 1, batch_index + 1))

# calculate the age of each individual assertion
def calculate_age(datestring):
    born = datetime.datetime.strptime(datestring, '%Y-%m-%d')
    today = datetime.date.today()
    return today.year - born.year - ((today.month, today.day) < (born.month,
                                                                 born.day))

# calculate age weight of each assertion
def age_weight(age):

    if age < 2:
        return 1
    elif 2 <= age <= 7:
        return ((12 - age) / 10)
    elif age > 7:
        return 0.5

# when using a list of ages, catch empty list and return average age
def average_list_age(vid, age_list):
    #numpy RuntimeWarning for empty age_list, suppress warning and log 
    with warnings.catch_warnings():
        warnings.simplefilter('error', category=RuntimeWarning)

        try:
            return np.nanmean(age_list)

        except:
            logging.warn('{} does not have valid clinical assertions!'
                          .format(vid))
            return None

# try a key to see if it is in the dict or warn
def key_test(test_dict, test_key):
    
    try:
        return test_dict[test_key]
    except KeyError:
        logging.warn('{} is not an expected key'.format(test_key))


class VariationClass:
    
    def __init__(self, variationreport):
        self.VID = variationreport.attrib['VariationID']
        self.CVVT = variationreport.attrib['VariationType']
        self.allele_parse(variationreport)
        self.observation_parse(variationreport)
        self.assertion_table_stats(variationreport)
        self.analysis_stats()

    # parse the Allele subtree of variation report
    def allele_parse(self, variationreport):
        RS = []
        Alt = []
        vcf_match = []
        
        for index, alleles in enumerate(variationreport.findall('./Allele')):
            
            try:
                RS.append(alleles.find('./XRefList/XRef[@DB="dbSNP"]')
                          .get('ID'))
            except:
                RS.append('')
            
            try:
                Alt.append(alleles
                           .find('./SequenceLocation[@Assembly="GRCh38"]')
                           .get('alternateAllele'))
            except:
                Alt.append('')
                
            vcf_match.append('{}|{}'.format(RS[index], Alt[index]))
            
        self.RSID = RS
        self.CVMA = Alt
        
        if len(vcf_match) > 1:
            self.vcfmatch = vcf_match
        else:
            self.vcfmatch = '.'

    # parse the PhenotypeList subtree of variation report
    def pheno_parse(self, observation):
        pheno_list = []
        
        for phenotype in observation.findall('./PhenotypeList/Phenotype'):
            pheno_list.append(phenotype.get('Name'))
        
        self.CVDS = ';'.join(pheno_list)

    # parse the ObservationList subtree of variation report
    def observation_parse(self, variationreport):
        run_already = False

        for observation in variationreport \
                .findall('./ObservationList/Observation'):
            
            if (observation.get('VariationID') == self.VID and
                    not run_already):
                run_already = True
                reviewstat = observation.find('ReviewStatus').text
                self.CVSZ = g.star_dict[reviewstat]
                self.CVCS = observation \
                    .find('./ClinicalSignificance/Description').text
                self.CVLE = observation.find('./ClinicalSignificance') \
                    .get('DateLastEvaluated')
                self.pheno_parse(observation)

            elif (observation.get('VariationID') == self.VID and
                    run_already):
                logging.warn('{} has multiple observation fields in its recor'
                             'd omitting as an annotation error. Check rsid(s'
                             ') {} manually'.format(self.VID, self.RSID))
                continue

            else:
                continue
        
    # parse the ClinicalAssertionList subtree of variation report
    def assertion_table_stats(self, variationreport):
        raw_score = []
        age_list = []

        for assertion in variationreport \
                .findall('./ClinicalAssertionList/GermlineList/Germline'):
            revstat_key = assertion.find('ReviewStatus').text
            score = key_test(g.cutoff, revstat_key)
            sigval_key = assertion.find('./ClinicalSignificance/Description') \
                .text
            sig_value = key_test(g.significance, sigval_key)
                
            if score > 0 and sig_value != 0:
                try:
                    age = calculate_age(assertion
                                        .find('./ClinicalSignificance')
                                        .get('DateLastEvaluated'))
                except:
                    logging.warn('{} has a missing assertion date!'
                                 .format(self.VID))
                    continue
                    
                age_list.append(age)
                raw_score.append(score * sig_value * age_weight(age))

        self.CTRS = sum(raw_score)
        self.CVNA = len(age_list)
        self.CTAA = average_list_age(self.VID, age_list)

        logging.debug('age list size: {}, raw_score size: {}'
                      .format(len(age_list), len(raw_score)))

    # calculating the analytical stats: CTWS, CTRR.
    def analysis_stats(self):
        sig_diff = {}
        
        for sig, value in g.ctws_cutoffs.items():
            sig_diff[sig] = abs(value - self.CTRS)
        
        self.CTWS = min(sig_diff, key=sig_diff.get)
        
        if self.CTWS in self.CVCS:
            self.CTRR = 'Consistent Classification'
        else:
            self.CTRR = 'Inconsistent result-under construction'
            
# test
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logging.debug('this is a module test')
    tree = ET.parse('../test/sample.xml')
    clinvarresult = tree.getroot()
    for var_index, variationreport in enumerate(clinvarresult):
        variant = VariationClass(variationreport)
        pprint.pprint(variant.__dict__)
    logging.debug('{} variant objects parsed in test'.format((var_index + 1)))
