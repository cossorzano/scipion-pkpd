"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-pkpd',  # Required
    version='1.0.0',  # Required
    description='scipion-pkpd.',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/cossorzano/scipion-pkpd',  # Optional
    author='Carlos Oscar Sorzano',  # Optional
    author_email='info@kinestat.com',  # Optional
    keywords='pkpd',  # Optional
    packages=find_packages(),
    package_data={  # Optional
       'pkpd': ['protocols.conf'],
    },
    install_requires=['scipy==0.15.0']
)
