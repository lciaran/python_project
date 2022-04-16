from setuptools import setup, find_packages

REQUIREMENTS = ['argparse', 'wget', 'gzip', 'shutil','biopython', 'statistics', 'numpy>=1.18.5', 'pandas', 'seaborn',  'matplotlib']

setup(name='ProtFLEXpreD',
      version='1.0',
      description='Program that calculates the flexibility of the proteins',
      author='Laura Ciaran Alfano and Neus Pou Amengual',
      author_email='laura.ciaran01@estudiant.upf.edu and neus.pou01@estudiant.upf.edu',
      url='https://github.com/lciaran/python_project',
      install_requires = REQUIREMENTS,
      packages = find_packages('ProtFLEXpreD'),
      package_dir = {"": "ProtFLEXpreD"},
      include_package_data= True,
      py_modules = ['ProtFLEXpreD', 'ProtFLEXpreD_functions', 'ProtFLEXpreD_graphical_representations', 'dictionaries'],
      python_requires='>=3',
      package_data={"":["*.md"]}
      )