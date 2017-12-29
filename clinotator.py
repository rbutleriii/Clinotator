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

__version__ = "0.0.1"

import argparse, re
import xml.etree.ElementTree as ET
import pandas as pd
import Bio.Entrez as Entz

# Args
parser = argparse.ArgumentParser(
    prog = 'Clinotator' ,
    description =
    '''
    Copyright (C) 2017  Robert R Butler III
    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
    '''
    )
parser.add_argument( "--version" , action = 'version' , version = '%(prog)s v' + __version__ )
parser.add_argument( "input" , metavar = ( '*.txt(vcf)' ) , help = "input file(s) (returns outfile for each)" , nargs = '+' )
parser.add_argument( "--type" , choices = [ 'vid' , 'rsid' , 'vcf' ] , nargs = 1 ,
                    help = '''
                    input file type: ClinVar Variation ID list
                                     dbSNP rsID list
                                     vcf file (output vcf generated)
                    ''' )
args = parser.parse_args()


