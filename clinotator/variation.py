#!/usr/bin/env python

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
import decimal
from math import fsum
from typing import List, Optional
import numpy as np
import xml.etree.ElementTree as ET
import global_vars as g


def query_parsing(
        variant_objects: List['VariationClass'],
        query_results: List[str]) -> None:
    """Parse each ClinVar XML batch into ``VariationClass`` objects.

    Args:
        variant_objects: List to be populated in-place with one
            ``VariationClass`` per variation report found across all
            batches.
        query_results: Raw ClinVar XML text, one string per batch, as
            produced by ``getncbi.get_ncbi_xml``.
    """
    for batch_index, batch in enumerate(query_results):
        clinvarresult = ET.fromstring(batch)
        for var_index, variationreport in enumerate(clinvarresult):
            variant = VariationClass(variationreport)
            variant_objects.append(variant)
        logging.debug('{} variant objects parsed in batch {}'
                      .format(var_index + 1, batch_index + 1))

def calculate_age(datestring: str) -> int:
    """Calculate the age in whole years since a given date.

    Args:
        datestring: A date in ``'%Y-%m-%d'`` format, e.g. the assertion's
            ``DateLastEvaluated``.

    Returns:
        The number of full years elapsed between ``datestring`` and today.
    """
    born = datetime.datetime.strptime(datestring, '%Y-%m-%d')
    today = datetime.date.today()
    return today.year - born.year - ((today.month, today.day) < (born.month,
                                                                 born.day))

def age_weight(age: int) -> float:
    """Compute the age-based weight applied to a clinical assertion.

    More recent assertions are weighted more heavily than older ones.

    Args:
        age: Age of the assertion in years, as from ``calculate_age``.

    Returns:
        ``1.0`` for age < 2, a linearly decreasing weight from just under
        1.0 down to 0.5 for age in [2, 7), and ``0.5`` for age >= 7.
    """
    if age < 2:
        return 1
    elif 2 <= age < 7:
        return ((11 - age) / 10)
    elif age >= 7:
        return 0.5

def average_list_age(vid: str, age_list: List[int]) -> Optional[float]:
    """Average a list of assertion ages, tolerating an empty list.

    Args:
        vid: The variant's ClinVar Variation ID, used only for logging.
        age_list: Ages (in years) of the variant's qualifying clinical
            assertions.

    Returns:
        The mean age, or ``None`` if ``age_list`` is empty.
    """
    #numpy RuntimeWarning for empty age_list, suppress warning and log 
    with warnings.catch_warnings():
        warnings.simplefilter('error', category=RuntimeWarning)

        try:
            return np.nanmean(age_list)

        except:
            logging.warning('VID: {} does not have valid clinical assertions!'
                          .format(vid))
            return None

def key_test(test_dict: dict, test_key: str):
    """Look up a key in a dict, logging and re-raising on ``KeyError``.

    Args:
        test_dict: The dict to look up ``test_key`` in (e.g.
            ``global_vars.cutoff`` or ``global_vars.significance``).
        test_key: The key to look up.

    Returns:
        The value at ``test_dict[test_key]``.

    Raises:
        KeyError: If ``test_key`` is not present in ``test_dict``.
    """
    try:
        return test_dict[test_key]
    except KeyError:
        logging.warning('{} is not an expected key'.format(test_key))
        raise KeyError

def reclassification_tree(ctps_index: int, cvcs_index: int) -> Optional[int]:
    """Score how strongly a variant's evidence suggests reclassification.

    Compares the bin index of the Clinotator-predicted significance
    (CTPS) against the bin index of ClinVar's submitted clinical
    significance (CVCS) to produce a reclassification-recommendation
    score.

    Args:
        ctps_index: Index of the variant's Clinotator Predicted
            Significance within ``global_vars.ctps_cutoffs``.
        cvcs_index: Index of the variant's ClinVar clinical significance
            within ``global_vars.ctps_cutoffs``.

    Returns:
        An integer score from 0 (no discordance) to 3 (highest-priority
        discordance), or ``None`` if the index difference and category
        combination don't match any defined case.
    """
    index_diff = abs(ctps_index - cvcs_index)
    
    if index_diff == 0:
        return 0
    elif index_diff == 1 and cvcs_index in (0,1,2,4,5,6) and ctps_index != 3:
        return 1
    elif index_diff == 1 and cvcs_index in (0,1,2,4,5,6) and ctps_index == 3:
        return 2
    elif index_diff == 1 and cvcs_index == 3:
        return 2
    elif (index_diff == 2 and cvcs_index in (0, 1, 2)
            and ctps_index in (0, 1, 2)):
        return 2
    elif (index_diff == 2 and cvcs_index in (0, 1, 2)
            and ctps_index not in (0, 1, 2)):
        return 3
    elif index_diff == 2 and cvcs_index == 3:
        return 3
    elif (index_diff == 2 and cvcs_index in (4, 5, 6)
            and ctps_index in (4, 5, 6)):
        return 2
    elif (index_diff == 2 and cvcs_index in (4, 5, 6)
            and ctps_index not in (4, 5, 6)):
        return 3
    elif index_diff > 2:
        return 3
                

class VariationClass:
    """Parsed ClinVar variation report with derived Clinotator metrics.

    Each instance represents a single ClinVar Variation ID, with
    attributes populated from the raw XML (e.g. ``VID``, ``CVCS``) as well
    as Clinotator's own derived scoring metrics (e.g. ``CTRS``, ``CTPS``).

    Attributes:
        VID (str): ClinVar Variation ID.
        CVVT (str): ClinVar variant type.
        CVSZ (int): ClinVar review-status star rating.
        rsID (List[str]): dbSNP rsID(s) associated with the variant.
        CVAL (List[str]): GRCh38 alternate allele(s).
        vcfmatch (Union[List[str], str]): ``'rsid|allele'`` pairs, or
            ``'.'`` if there is only one allele.
        CVCS (str): ClinVar clinical significance description.
        CVLE (str): Date the ClinVar interpretation was last evaluated,
            or ``'.'`` if not available.
        CVDS (str): Semicolon-joined disease/phenotype terms with their
            per-assertion significance codes.
        CVNA (int): Number of qualifying (germline, criteria-provided)
            clinical assertions.
        CTAA (Optional[float]): Average age (years) of qualifying clinical
            assertions.
        CTRS (Optional[float]): Clinotator raw weighted pathogenicity
            score, or ``None`` if fewer than 2 qualifying assertions.
        CTPS (Optional[str]): Clinotator-predicted clinical significance
            bin.
        CTRR (Union[int, str]): Reclassification recommendation score, or
            ``'.'`` if not computable.
    """

    def __init__(self, variationreport: ET.Element) -> None:
        """Parse one ClinVar ``VariationArchive`` element into attributes.

        Args:
            variationreport: An ``xml.etree.ElementTree.Element`` for a
                single ClinVar variation report (one child of the batch
                root element).
        """
        self.VID = variationreport.attrib['VariationID']
        self.CVVT = variationreport.attrib['VariationType']
        revstat = variationreport.find('.InterpretedRecord/ReviewStatus').text
        self.CVSZ = g.star_dict[revstat]
        self.allele_parse(variationreport)
        self.observation_parse(variationreport)
        self.assertion_table_stats(variationreport)
        self.analysis_stats()

    def haplos(self, variationreport: ET.Element) -> str:
        """Choose the XPath for this report's alleles, handling haplotypes.

        Args:
            variationreport: The variation report element for this
                variant.

        Returns:
            An XPath string locating the ``SimpleAllele`` element(s),
            accounting for the extra ``Haplotype`` level present when
            ``self.CVVT == "Haplotype"``.
        """
        if self.CVVT == "Haplotype":
            all_term = '.InterpretedRecord/Haplotype/SimpleAllele'
        else:
            all_term = '.InterpretedRecord/SimpleAllele'
        return all_term

    def allele_parse(self, variationreport: ET.Element) -> None:
        """Parse allele-level fields (rsID, ALT allele) from the report.

        Populates ``self.rsID``, ``self.CVAL``, and ``self.vcfmatch``.

        Args:
            variationreport: The variation report element for this
                variant.
        """
        RS = []
        Alt = []
        vcf_match = []
        
        for index, alleles in enumerate(variationreport
                .findall(self.haplos(variationreport))):
            
            try:
                RS.append(alleles.find('./XRefList/XRef[@DB="dbSNP"]')
                        .get('ID'))
            except:
                RS.append('.')
            
            try:
                Alt.append(alleles
                        .find('./Location/SequenceLocation[@Assembly="GRCh38"]')
                        .get('alternateAlleleVCF'))
            except:
                Alt.append('.')
                
            vcf_match.append('{}|{}'.format(RS[index], Alt[index]))
            
        logging.debug('RS list -> {}\nAlt list -> {}'.format(RS, Alt))
        self.rsID = RS
        self.CVAL = Alt
        
        if len(vcf_match) > 1:
            self.vcfmatch = vcf_match
        else:
            self.vcfmatch = '.'

    def observation_parse(self, variationreport: ET.Element) -> None:
        """Parse the clinical-significance interpretation for the variant.

        Populates ``self.CVCS`` and ``self.CVLE``. Warns and skips
        additional clinical-significance interpretation blocks if more
        than one is present, since that indicates an annotation error.

        Args:
            variationreport: The variation report element for this
                variant.
        """
        run_already = False

        for interpretation in variationreport \
                .findall('./InterpretedRecord/Interpretations/Interpretation'):
            
            if (interpretation.get('Type') == "Clinical significance" and
                    not run_already):
                run_already = True
                self.CVCS = interpretation \
                    .find('./Description').text
                try:
                    self.CVLE = interpretation.attrib['DateLastEvaluated']
                except KeyError as e:
                    self.CVLE = '.'
                    logging.warning('VID {} doesn\'t have a DateLastEvaluated!'
                                    .format(self.VID))

            elif interpretation.get('VariationID') == self.VID and run_already:
                logging.warning('{} has multiple interpretation fields in its '
                                'record omitting as an annotation error. Check'
                                ' rsid(s) {} manually'.format(self.VID,
                                                              self.rsID))
                continue

            else:
                continue
        
    def pheno_parse(self, assertion: ET.Element, sig_key: str) -> List[str]:
        """Build disease/phenotype labels for one clinical assertion.

        Args:
            assertion: A ``ClinicalAssertion`` element from the variation
                report.
        sig_key: The short significance code (e.g. ``'P'``, ``'LB'``) to
                annotate each phenotype term with.

        Returns:
            A list of ``'DB:ID(sig_key)'`` strings, one per associated
            phenotype/trait; or a single ``'Not_Provided(sig_key)'`` entry
            if the assertion lists no phenotype cross-references.
        """
        pheno_list = []
        
        for phenotype in assertion.findall('./TraitSet/Trait/XRef'):
            d_name = ':'.join([phenotype.attrib['DB'], phenotype.attrib['ID']])
            pheno_list.append('{}({})'.format(d_name, sig_key))
        
        if not pheno_list:
            pheno_list.append('{}({})'.format("Not_Provided", sig_key))
        logging.debug('Disease list for {}: {}'.format(assertion.attrib['ID'],
                                                       pheno_list))
        return pheno_list

    def assertion_table_stats(self, variationreport: ET.Element) -> None:
        """Score qualifying clinical assertions and aggregate the results.

        Iterates over each ``ClinicalAssertion`` for this variant,
        filters to germline-origin assertions with valid review status
        and significance, and computes a weighted raw score per
        assertion. Populates ``self.CVDS``, ``self.CVNA``, ``self.CTAA``,
        and ``self.CTRS``.

        Args:
            variationreport: The variation report element for this
                variant.
        """
        raw_score = []
        age_list = []
        cvds_list = []
        
        for assertion in variationreport.findall(
                './InterpretedRecord/ClinicalAssertionList/ClinicalAssertion'):
            observ_set = {"germline", "de novo", "maternal", "paternal",
                          "inherited", "unknown", "uniparental", "biparental"}
            observ_list = {x.text.lower() for x in assertion
                    .findall('./ObservedInList/ObservedIn/Sample/Origin')}
            logging.debug('Origin List for {}: {}'
                          .format(assertion.attrib['ID'], observ_list))
            try:
                assert len(observ_set.intersection(observ_list)) > 0
                revstat_key = assertion.find('ReviewStatus').text
                score = key_test(g.cutoff, revstat_key)
                sigval_key = assertion.find('./Interpretation/Description') \
                    .text
                try:
                    sig_value = key_test(g.significance, sigval_key)
                except:
                    logging.warn('Assertion {} for VID {} is incorrectly forma'
                            'tted'.format(assertion.attrib['ID'], self.VID))
                    continue
    
                if score > 0 and sig_value[0] != 0:
                    try:
                        age = calculate_age(assertion.find('./Interpretation')
                                            .get('DateLastEvaluated'))
                    except:
                        logging.debug('Assertion {} for VID {} is missing an a'
                                      'ssertion date!'.format(
                                      assertion.attrib['ID'],self.VID))
                        continue
                        
                    age_list.append(age)
                    D = decimal.Decimal
                    raw_score.append(float(D(str(score)) * D(str(sig_value[0]))
                                     * D(str(age_weight(age)))))
                    logging.debug('score: {} sig_value: {} age_weight: {} age:'
                                  ' {}'.format(score, sig_value[0],
                                               age_weight(age), age))
    
                    cvds_list += self.pheno_parse(assertion, sig_value[1])
            except AssertionError as a:
                logging.debug('no germline reports for assertion {}, skipping'
                              .format(assertion.attrib['ID']))
                continue

        self.CVDS = ';'.join(cvds_list)
        if not self.CVDS:
            self.CVDS = '.'
        
        self.CVNA = len(age_list)
        self.CTAA = average_list_age(self.VID, age_list)
        # logging.debug('VID: {} -> age list size: {}, raw_score size: {}'
        #               .format(self.VID, len(age_list), len(raw_score)))

        if len(raw_score) >= 2:
            # logging.debug('raw_score list {}'.format(raw_score))
            self.CTRS = fsum(raw_score)
        else:
            self.CTRS = None
            logging.debug('VID: {} has fewer than 2 complete assertions, no C'
                          'TRS will be generated'.format(self.VID))

    def analysis_stats(self) -> None:
        """Compute the Clinotator Predicted Significance and CTRR.

        Requires ``self.CTRS``/``self.CVNA`` to already be set (by
        ``assertion_table_stats``). Populates ``self.CTPS`` and
        ``self.CTRR``, or sets both to their "not computable" values
        (``None`` and ``'.'`` respectively) if there are too few
        assertions or the ClinVar significance doesn't map to a known
        bin.
        """
        if self.CVNA < 2:
            self.CTPS = None
            self.CTRR = '.'
            return
        
        # logging.debug('CTRS score: {}'.format(self.CTRS))
        bins = [x[1] for x in g.ctps_cutoffs[:6]]
        bins = np.nextafter(bins, np.minimum(bins,0))
        ctps_index = np.digitize(self.CTRS, bins, right=True).tolist()
        self.CTPS = g.ctps_cutoffs[ctps_index][0]
        
        run_already = False
        cvcs_index = None
        for clinsig in self.CVCS.split(', '):
            if clinsig in [x[0] for x in g.ctps_cutoffs] and not run_already:
                run_already = True
                cvcs_index = [x[0] for x in g.ctps_cutoffs].index(clinsig)
            elif clinsig == 'Conflicting interpretations of pathogenicity':
                run_already = True
                cvcs_index = 3
            elif clinsig in [x[0] for x in g.ctps_cutoffs] and run_already:
                logging.warning('multiple CVCS statements for {}, only using '
                             'first one!'.format(self.VID))
        
        if cvcs_index is None:
            logging.warning('ClinVar significance for {} does not include B,B/'
                            'LB,LB,US,LP,LP/P,P'.format(self.VID))
            self.CTPS = None
            self.CTRR = '.'
            return
        
        self.CTRR = reclassification_tree(ctps_index, cvcs_index)
        
# test
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logging.debug('this is a module test')
    tree = ET.parse('../test/sample2.xml')
    clinvarresult = tree.getroot()
    for var_index, variationreport in enumerate(clinvarresult):
        variant = VariationClass(variationreport)
        pprint.pprint(variant.__dict__)
    logging.debug('{} variant objects parsed in test'.format((var_index + 1)))
