from setuptools import setup, find_packages
import os

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open(os.path.join('./', 'VERSION')) as version_file:
    version = version_file.read().strip()

setup(
    name='clinotator',
    version=version,
    description='Clinical interpretation of ambiguous ClinVar annotations',
    long_description=readme,
    author='Robert R Butler III',
    author_email='rbutler@northshore.org',
    url='https://github.com/rbutleriii/clinotator',
    license=license,
    packages=find_packages(),
    install_requires=['numpy>=1.8.0', 'pandas>=0.20.0', 'biopython>=1.70']
)

