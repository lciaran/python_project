from setuptools import setup

REQUIREMENTS = ['argparse', 'wget','biopython', 'statistics', 'numpy>=1.18.5', 'pandas', 'seaborn',  'matplotlib']

setup(name='ProtFLEXpreD',
      version='1.0',
      description='Program that calculates the flexibility of the proteins',
      author='Laura Ciaran Alfano and Neus Pou Amengual',
      author_email='laura.ciaran01@estudiant.upf.edu and neus.pou01@estudiant.upf.edu',
      url='https://github.com/lciaran/python_project',
      install_requires = REQUIREMENTS,
      packages = ['ProtFLEXpreD', 'ProtFLEXpreD'],
      py_modules = ['ProtFLEXpreD', 'PFD_functions', 'PFD_representations', 'dictionaries'],
      python_requires='>=3.8'
      )