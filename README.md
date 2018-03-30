[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1210203.svg)](https://doi.org/10.5281/zenodo.1210203)

# Clinotator
## Synopsis

### Clinical interpretation of ambiguous ClinVar annotations

This project takes variants as input and queries NCBI eutilities to generate ClinVar Variation Report<sup>1</sup> scoring metrics. The overall goal is to generate annotations of use for given batches of variants to inform clinical interpretation. The metrics include:

*	Clinotator Raw Score -  A weighted metric of pathogenicity based on submitter type, assertion type and assertion age. 
*	Average Clinical Assertion Age -  Clinical assertions with criteria provided are counted, and their average age is calculated.
*	Clinotator Predicted Significance -  A predicted clinical significance based on the weighted distribution of all two star variants in ClinVar with two or more clinical assertions.
*	Reclassification Recommendation -  A ranking of the impact of reclassification based on the Clinotator Predicted Significance. 

These and other stats are returned on a per variant basis, in a table, and additionally as vcf annotations if a vcf file is provided. See below for more information. 

## Code Example

```bash
usage: clinotator.py [-h] [--log] [--long-log] [-o prefix] [--version] -e
                     EMAIL -t {vid,rsid,vcf}
                     file [file ...]

Clinical interpretation of ambiguous ClinVar annotations

positional arguments:
  file               input file(s) (returns outfile for each)

optional arguments:
  -h, --help         show this help message and exit
  --log              create logfile
  --long-log         create detailed logfile
  -o prefix          choose an alternate prefix for outfiles
  --version          show program's version number and exit

required arguments:
  -e EMAIL           NCBI requires an email for database queries
  -t {vid,rsid,vcf}  vid - ClinVar Variation ID list
                     rsid - dbSNP rsID list
                     vcf - vcf file (output vcf generated)
```

### Minimum Inputs

Three required bits of information: **(1)** the type of input file, **(2)** the file itself and **(3)** your email address. The input file should be either a single column of IDs in a text file, or a vcf. If a vcf(s) is selected, clinotator will generate an annotated output file. In all cases a tab-delimited table file will be produced. The email is required by NCBI/biopython as NCBI enforces a strict 3 queries per second limit. Before they ban your IP address they will attempt to contact you, but that should not be a problem as clinotator uses batch queries and the history server.

### Optional Arguments

Additional arguments include a log file (--log), a more detailed log file (--long-log) and specification of the output file prefix (-o, the default is clinotator). Also help and version messages are available.

## Motivation

While ClinVar has become an indispensable resource for clinical variant interpretation, its sophisticated structure provides a daunting learning curve. Often the sheer depth of types of information provided can make it difficult to analyze variant information with high throughput. Clinotator is a fast and lightweight tool to extract important aspects of criteria-based clinical assertions and uses that information to generate several metrics to assess the strength and consistency of the evidence supporting the variant clinical significance. Clinical assertions are weighted by significance type, age of submission and submitter expertise category to filter outdated or incomplete assertions that otherwise confound interpretation. This can be accomplished in batches: either lists of Variation IDs or dbSNP rsIDs, or with vcf files which are additionally annotated. Clinotator slices out problem variants in minutes without extensive computational effort—just using a personal computer. With the rapidly growing body of variant evidence, most submitters and researchers have limited resources to devote to variant curation. Clinotator provides efficient, systematic prioritization of discordant variants in need of reclassification. The hope is that this tool can inform ClinVar curation and encourage submitters to keep their clinical assertions current by focusing their efforts. Additionally, researchers can utilize new metrics to analyze variants of interest in pursuit of new insights into pathogenicity.

## Installation

Implemented in python (tested on 2.7.12 and >=3.4). You can `git clone` or download the zipfile and unpack. Add the folder location to your ~/.bash_profile or:

```
export PATH=$PATH:path/to/folder/Clinotator/clinotator
``` 

Examples of each input file type are provided in the test subfolder. For instance:

```
cd path/to/Clinotator/test
clinotator.py -t vid -e A.N.Other@example.com test.vid
```

Should produce the following warnings and a clinotator.test.tsv file:

```
INFO:root:Starting on test.vid
INFO:root:Going to download record 1 to 13
INFO:root:Download time: 0.026888549999997242 min, Batches run -> 1
WARNING:root:128294 has a missing assertion date!
WARNING:root:128297 has a missing assertion date!
WARNING:root:ClinVar significance for 3521 does not include B,B/LB,LB,US,LP,LP/P,P
WARNING:root:VID: 55794 does not have valid clinical assertions!
```

The warnings, as well as some additional information can be stored in the log file with `--log`. `--long-log` will store detailed debugging information, but the file will be larger than the output tsv file. Both log files append information, so batch runs or especially large lists of variants may result in large file sizes.

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

### Memory/System requirements

Clinotator was designed in a Linux environment and implemented in Python (2.7 or >=3.4), and can run in similar OSX and possibly Windows Python environments. It can be run on a personal computer with relatively modest system requirements; a minimum of 2GB available RAM.  

As Clinotator keeps the NCBI xml results in memory, there is a significant memory usage. At the time of writing, the entire ClinVar xml set is approaching 6GB. Loading the entire set into memory is doable with at least 8GB of memory, though it is recommended that you batch your queries in this rare case. More typical usage for subsets of ClinVar or batch vcf annotations should not pose a memory issue.

## Details on metrics

### ClinVar Metrics

<dl>
	<dt>ClinVar Clinical Significance (CVCS)</dt>
	<dd>Clinical significance reported by ClinVar.<sup>2</sup> Ratings metrics are based on the five ACMG/AMP recommended <sup>3</sup> classifications for Mendelian disorders: Benign, Likely benign, Uncertain significance, Likely pathogenic and Pathogenic. Other Clinical significance values are reported, but not factored into the Clinotator metrics.</dd>
</dl>
<dl>
	<dt>ClinVar Stars (CVSZ)</dt>
	<dd>Star rating given by clinvar. Ranges from zero to four.<sup>4</sup></dd>
</dl>
<dl>
	<dt>ClinVar Number of Clinical Assertions (CVNA)</dt>
	<dd>The number of Clinvar Submissions possessing a clinical assertion (with criteria provided). This measure excludes submissions without assertion criteria, including "literature reviews", which are a type of evidence as opposed to an assertion. Additionally, submitter assertions without defined criteria are also omitted. Most assertions with criteria meet or exceed the guidelines put forth by the American College of Medical Genetics and Genomics (ACMG) and the Association for Molecular Pathology (AMP) in 2015.<sup>3</sup></dd>
</dl>
<dl>
	<dt>ClinVar Conditions/Diseases (CVDS)</dt>
	<dd>Conditions reported to be associated with this variant.</dd>
</dl>
<dl>
	<dt>ClinVar Last Evaluated (CVLE)</dt>
	<dd>The date the clinical significance of the variation report was last evaluated. Note this is not the date the variation report was last updated, but the date in the <ClinicalAssertionList> field of the ClinVar xml connected to the Review Status.</dd>
</dl>
<dl>
	<dt>ClinVar Variant Type (CVVT)</dt>
	<dd>The type of variation in ClinVar. Currently defined as either "Simple" with a single AlleleID or "Haplotype" if multiple AlleleIDs are involved.</dd>
</dl>

### Clinotator Metrics

<dl>
	<dt>Clinotator Raw Score (CTRS)</dt>
	<dd>A weighted metric of pathogenicity based on submitter type, assertion type and assertion age. The type of submitter is weighted based on expertise, with regular clinical assertions unweighted at 1.00, expert reviewers receiving a 1.10 and practice guidelines receiving a score of 1.25.<br /><br /> The age of the assertion is weighted as new data is incorporated into assertions as well as previous data, creating a larger set of evidence over time. For the first two years, there is no weight, then there is a 10% reduction in weight per year through 6 years , at which point the penalty stays at a static 50% weight thereafter.<br /><br /> The assertion type is that largest weight, with values of: Benign(B) = -6, Likely benign(LB) = -3, Uncertain significance(US) = -0.3, Likely pathogenic(LP) = 3 and Pathogenic(P) = 6. For more information on the weighting decisions, see our publication.<sup>5</sup></dd>
</dl>
<dl>
	<dt>Average Clinical Assertion Age (CTAA)</dt>
	<dd>As described above, the clinical assertions with criteria provided are counted, and their average age is calculated.</dd>
</dl>
<dl>
	<dt>Clinotator Predicted Significance (CTPS)</dt>
	<dd>This is a *predicted* clinical significance based on the weighted distribution of all variants in ClinVar with two or more clinical assertions (as of a Clinotator version release date). The ratings are calculated as previously described, on nonparametric prediction intervals with a given confidence of classification. See Figure 1 in our publication for details.<sup>5</sup></dd>
</dl>
<dl>
	<dt>Clinotator Reclassification Recommendation (CTRR)</dt>
	<dd>This field ranks reclassification priority based on the difference between the CVCS and the CTWS. This field only includes the seven values of ClinVar clinical significance associated with Mendelian diseases (B, B/LB, LB, US/CI, LP, LP/P, P). For the purposes of reclassification, "Conflicting interpretations of pathogenicity" is scored the same as Uncertain significance.</dd>
</dl>

*	0 - Reclassification unlikely, consistent identity or insufficient information for a recommendation
*	1 - Low priority reclassification, minor change without clinical impact
*	2 - Medium priority reclassification, minor change of some clinical impact
*	3 - High priority reclassification, significant change in clinical impact

<dl>
	<dd>For a detailed decision schema, see our publication (see Fig 2).<sup>5</sup></dd>
</dl>

## Citation

Citation

'''
BibTex format citation coming soon
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
<sup>3</sup>  Richards, S., et al. (2015). "Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology." Genet Med 17(5): 405-424.  
<sup>4</sup>  https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/  
<sup>5</sup>  Our paper, Coming Soon  
