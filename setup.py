from setuptools import setup, find_packages
import clinotator

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='clinotator',
    version=clinotator.__version__,
    description='Clinical interpretation of ambiguous ClinVar annotations',
    long_description=readme,
    author='Robert R Butler III',
    author_email='rbutler@northshore.org',
    url='https://github.com/rbutleriii/clinotator',
    license=license,
    packages=find_packages(exclude=())
)

