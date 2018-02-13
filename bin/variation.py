#!/usr/bin/env python3

'''
Below are the functions for the Variation file object which contains variant
info for metrics.

Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''


import logging
import xml.etree.ElementTree as ET


class VariationClass:
    
    def __init__(self, variationreport):
        self.variationreport = variationreport #this is a ET Elemenmt object
    
    # # # # takes python object from clinvar and puts info into dataframe
    def parse_to_tbl(self):
        parse_submissions(submission_subtree)
    
    # # # # takes subtree from main for submitters and calculates submission stats
    def parse_submissions(self):
        # build dataframe
        # aggregate functions for data
        pass
    
    def clinvar_stats(self):
        self.var_id
        self.rsid
        self.clinsig
        self.starz
        self.num_assert
        self.disease
        self.last_update
   
    def ct_score():
        pass
    
    def ct_average_age():
        pass
    
    def ct_weighted_significance():
        pass
    
    def ct_reclassification():
        pass
