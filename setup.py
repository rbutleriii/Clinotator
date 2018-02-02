from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='clinotator',
    version='0.1.0',
    description='add desc later',
    long_description=readme,
    author='Robert R Butler III',
    author_email='rbutler@northshore.org'
    url='https://github.com/rbutleriii/clinotator',
    license=license,
    packages=find_packages(exclude=())
)

