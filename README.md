# Clinotator
## Synopsis

### Clinical interpretation of ambiguous ClinVar annotations

This project takes variants as input and queries NCBI eutilities to generate scoring metrics. The overall goal is to generate annotations of use for given batches of variants to inform clinical interpretation. The metrics include:

*	Clinotator Score - ...
*	Average Clinical Assertion Age - ...
*	Clinotator Weighted Significance - ...
*	Reclassification Recommendation - ...

These and other stats are returned on a per variant basis, in a table, and additionally as vcf annotations if a vcf file is provided. See below for more information. 

## Code Example

```css
usage: clinotator [-h] [--version] -t {vid,rsid,vcf} -e EMAIL file [file ...]

Clinical interpretation of ambiguous ClinVar annotations

positional arguments:
  file               input file(s) (returns outfile for each)

optional arguments:
  -h, --help         show this help message and exit
  --version          show program's version number and exit

required arguments:
  -t {vid,rsid,vcf}  vid - ClinVar Variation ID list
                     rsid - dbSNP rsID list
                     vcf - vcf file (output vcf generated)
  -e EMAIL           NCBI requires an email for querying their databases
```

### Minimum Inputs

Three required bits of information: the type of input file, the file itself and your email address. The input comes as either a simple text list file of IDs, or a vcf. If a vcf(s) is selected, clinotator will generate an annotated output file. In all cases a tab-delimited table file will be produced.


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
pipenv graph
biopython==1.70
  - numpy [required: Any, installed: 1.14.0]
pandas==0.22.0
  - numpy [required: >=1.9.0, installed: 1.14.0]
  - python-dateutil [required: >=2, installed: 2.6.1]
    - six [required: >=1.5, installed: 1.11.0]
  - pytz [required: >=2011k, installed: 2018.3]
```

The python setup.py will attempt to use setuptools to install dependencies. Likely choose your own use case.

### Proxy

biopython uses urllib to query eutilities (I believe now https only), so if you are behind a proxy add the following to your ~/.bash_profile:

```
export HTTP_PROXY=http://username:password@proxy.mydomain.com:8080
export HTTPS_PROXY=http://username:password@proxy.mydomain.com:8080
```

## Details on metrics

### ClinVar Metrics

ClinVar Clinical Significance (**CTCS**)
: ...
ClinVar Stars (**CTSZ**)
: ...
ClinVar \# of Clinical Assertions (**CTNA**)
: ...
ClinVar Diseases (**CTDS**)
: ...

### Clinotator Metrics

Clinotator Raw Score (**CTRS**)
: ...
Average Clinical Assertion Age (**CTAA**)
: ...
Clinotator Weighted Significance (**CTWS**)
: ...
Reclassification Recommendation (**CTRR**)
: ...

## Tests

For the main package, the following tests check the installation:

```
```

## License

Copyright (C) 2017  Robert R Butler III

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

## References

