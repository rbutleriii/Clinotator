#!/usr/bin/env python3

'''
Clinotator - Clinical interpretation of ambiguous ClinVar annotations
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
'''

__version__ = "0.1.0"

import argparse, re
import xml.etree.ElementTree as ET
import pandas as pd
import Bio.Entrez as Entz
Entz.tool = 'Clinotator' # preferred by NCBI
Entz.email = args.email

# Args
parser = argparse.ArgumentParser(
    prog = 'Clinotator' ,
    formatter_class = argparse.RawTextHelpFormatter ,
    description =
    '''
    Copyright (C) 2017  Robert R Butler III
    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>\
    ''')
parser.add_argument( "--version" , action = 'version' , version = '%(prog)s v' + __version__ )
parser.add_argument( '--email', nargs = 1, help = 'NCBI requires an email for querying their databases')
parser.add_argument( "input" , metavar = ( '*.file' ) , help = "input file(s) (returns outfile for each)" , nargs = '+' )
parser.add_argument( "--type" , choices = [ 'vid' , 'rsid' , 'vcf' ] , nargs = 1 ,
                    help =
'''\
input file type: vid - ClinVar Variation ID list
                 rsid - dbSNP rsID list
                 vcf - vcf file (output vcf generated)\
''')
args = parser.parse_args()

# Functions
def rsid_to_vid(rsid_list):
    # import Bio.Entrez as Entz
    # Entz.email = args.email
    # Entz.tool = 'Clinotator' # preferred by NCBI

    results = Entz.read(Entz.epost('snp' , id = ','.join(rsid_list)))
    webenv = results['WebEnv']
    query_key = results['QueryKey']
    
    return
    


