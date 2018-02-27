# Clinotator
## Synopsis

### Clinical interpretation of ambiguous ClinVar annotations

This project takes variants as input and queries NCBI eutilities to generate ClinVar Variation Report[1] scoring metrics. The overall goal is to generate annotations of use for given batches of variants to inform clinical interpretation. The metrics include:

*	Clinotator Raw Score -  A weighted metric of pathogenicity based on submitter type, assertion type and assertion age. 
*	Average Clinical Assertion Age -  The average age of clinical assertions made about a variant.
*	Clinotator Weighted Significance -  A predicted clinical significance based on prediction intervals of the Clinotator Raw Score.
*	Reclassification Recommendation -  A ranking of the impact of reclassification based on the Clinotator Weighted Significance. 

These and other stats are returned on a per variant basis, in a table, and additionally as vcf annotations if a vcf file is provided. See below for more information. 

## Code Example

```bash
usage: clinotator.py [-h] [--log] [-o prefix] [--version] -e EMAIL -t {vid,rsid,vcf} file [file ...]

Clinical interpretation of ambiguous ClinVar annotations

positional arguments:
  file               input file(s) (returns outfile for each)

optional arguments:
  -h, --help         show this help message and exit
  --log              create logfile
  -o prefix          choose an alternate prefix for outfiles
  --version          show program's version number and exit

required arguments:
  -e EMAIL           NCBI requires an email for querying their databases
  -t {vid,rsid,vcf}  vid - ClinVar Variation ID list
                     rsid - dbSNP rsID list
                     vcf - vcf file (output vcf generated)
```

### Minimum Inputs

Three required bits of information: (1) the type of input file, (2) the file itself and (3) your email address. The input comes as either a text ID list file, or a vcf. If a vcf(s) is selected, clinotator will generate an annotated output file. In all cases a tab-delimited table file will be produced. The email is required by NCBI/biopython as NCBI enforces a strict 3 queries per second limit. Before they ban your IP address they will attempt to contact you, but that should not be a problem as clinotator uses batch queries and the history server.

### Optional Arguments

Additional arguments include a log file and specification of the output file prefix.

## Motivation

While ClinVar has become an indispensable resource for clinical variant interpretation, that can be a double-edged sword. Often the sheer wealth of types of information can make it difficult to make a final interpretation. The sophisticated architecture the records also makes programmatic analysis of batches of variants challenging. This software filters information types for each variant to focus on clinical assertions being made about the variant. It generates several metrics by which to gauge the robustness of the overall clinical assertion and measure the variance in Interpretation. This can be done by batches of Variation IDs, batches of dbSNP rsIDs or by analysis and annotation of .vcf files. The hope is that this will help identify variants that are candidates for reclassification, and prioritize variants for further research.

## Installation

Implemented in python3 (tested >=3.5). You can `git clone` or download the zipfile and unpack as you like. Add the location to your ~/.bash_profile or:

```
export PATH=$PATH:path/to/folder/clinotator
``` 

### Dependencies

The following from my pipenv:

```
biopython==1.70
  - numpy [required: Any, installed: 1.14.0]
pandas==0.22.0
  - numpy [required: >=1.9.0, installed: 1.14.0]
  - python-dateutil [required: >=2, installed: 2.6.1]
    - six [required: >=1.5, installed: 1.11.0]
  - pytz [required: >=2011k, installed: 2018.3]
```

Numpy *should* work >= 1.9.0 and pandas >= 0.20.0, but install more recent versions if possible.

## Details on metrics

### ClinVar Metrics

<dl>
	<dt>ClinVar Clinical Significance (**CVCS**)</dt>
	<dd>Clinical significance reported by ClinVar.<sup>2</sup> Ratings metrics are based on the five ACMG/AMP recommended classifications for Mendelian disorders: Benign, Likely benign, Uncertain significance, Likely pathogenic and Pathogenic. Other Clinical significance values are reported, but not factored into the Clinotator metrics.</dd>
</dl>
<dl>
	<dt>ClinVar Stars (**CVSZ**)</dt>
	<dd>Star rating given by clinvar. Ranges from zero to four.<sup>3</sup></dd>
</dl>
<dl>
	<dt>ClinVar Number of Clinical Assertions (**CVNA**)</dt>
	<dd>The number of Clinvar Submissions possessing a clinical assertion (with criteria provided). This measure excludes submissions without assertion criteria, including "literature reviews", which are a type of evidence as opposed to an assertion. Additionally, submitter assertions without defined criteria are also omitted. Most assertions with criteria meet or exceed the guidelines put for by the American College of Medical Genetics and Genomics (ACMG) in 2013 and amended in 2015.<sup>4,5</sup></dd>
</dl>
<dl>
	<dt>ClinVar Conditions/Diseases (**CVDS**)</dt>
	<dd>Conditions reported to be associated with this variant.</dd>
</dl>
<dl>
	<dt>ClinVar Last Evaluated (**CVLE**)</dt>
	<dd>The date the clinical significance of the variation report was last evaluated. Note this is not the date the variation report was last updated, but the date in the \<ClinicalAssertionList\> field of the ClinVar xml connected to the Review Status.</dd>
</dl>
<dl>
	<dt>ClinVar Variant Type (**CVVT**)</dt>
	<dd>The type of variation in ClinVar. Currently defined as either "Simple" with a single AlleleID or "Haplotype" if multiple AlleleIDs are involved.</dd>
</dl>

### Clinotator Metrics

<dl>
	<dt>Clinotator Raw Score (**CTRS**)</dt>
	<dd>A weighted metric of pathogenicity based on submitter type, assertion type and assertion age. The type of submitter is weighted based on expertise, with regular clinical assertions unweighted at 1.00, expert reviewers receiving a 1.10 and practice guidelines receiving a score of 1.25.  

The age of the assertion is penalized as new data is incorporated into newer assertions as well as previous data, creating a larger set of evidence over time. For the first two years, there is no penalty, then there is a 10% reduction gradation in weight per year through 6 years , at which point the penalty stays at a static 50% reduction thereafter.  

The assertion type is that largest weight, with values of: Benign(B) = -5, Likely benign(LB) = -3, Uncertain significance(US) = -0.3, Likely pathogenic(LP) = 1.6 and Pathogenic(P) = 2.9. For more information on the weighting decisions, see our publication.<sup>6</sup></dd>
</dl>
<dl>
	<dt>Average Clinical Assertion Age (**CTAA**)</dt>
	<dd>As described above, the clinical assertions with criteria provided are counted, and their average age is calculated.</dd>
</dl>
<dl>
	<dt>Clinotator Predicted Significance (**CTWS**)</dt>
	<dd>This is a *predicted* clinical significance based on the weighted distribution of all variants in ClinVar with two or more clinical assertions (as of a Clinotator version release date). The ratings are calculated as previously described, on nonparametric prediction intervals with a 99% confidence of the given classification. See Figure 1 in our publication for details.<sup>6</sup></dd>
</dl>
<dl>
	<dt>Clinotator Reclassification Recommendation (**CTRR**)</dt>
	<dd>This field ranks reclassification priority based on the difference between the CVCS and the CTWS. This field only includes the seven values of clinical significance associated with Mendelian diseases (B, B/LB, LB, US/CI, LP, LP/P, P). For the purposes of reclassification, "Conflicting interpretations of pathogenicity" is scored the same as Uncertain significance.</dd>
</dl>

*	0 - Reclassification unlikely, consistent identity or 
*	1 - Low priority reclassification, change of low impact
*	2 - Medium priority reclassification, minor change of clinical impact
*	3 - High priority reclassification, significant change in clinical impact

<dl>
	<dd>For a detailed decision tree, see Figure 2 in our publication.<sup>6</sup></dd>
</dl>

## Citation

Citation

'''
BibTex format citation
'''


## License

Copyright (C) 2017  Robert R Butler III

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License v3 as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

## References

<sup>1</sup>  https://www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/  
<sup>2</sup>  https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/  
<sup>3</sup>  https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/  
<sup>4</sup>  ACMG Guidelines 2013  
<sup>5</sup>  ACMG Guidelines 2015  
<sup>6</sup>  Our paper  

