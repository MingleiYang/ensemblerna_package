import sys
from setuptools import setup, find_packages

setup(
      name='ensemblerna',
      version='1.0',
      author='Chanin Tolson Woods',
      author_email='laederachlab@gmail.com',
      description='A package for the visualization and comparison of RNA structural ensembles',
      license='GPL (>= 2)',
      keywords='visualization rna structural ensemble',
      url='https://www.ribosnitch-ensemblerna.rhcloud.com',
      packages=find_packages(exclude=['examples']),
      install_requires=['numpy', 'matplotlib', 'scipy', 'scikit-learn', 'jinja2', 'IPython', 'mpld3'],
      scripts=['bin/ensemblerna'],
      zip_safe=False
)