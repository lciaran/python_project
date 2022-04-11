from distutils.core import setup
PACKAGES = []
REQUIREMENTS = ['sys', 'os', 'argparse', 're', 'wget', 'biopython', 'statistics', 'numpy', 'pandas', 'seaborn',  'matplotlib']

setup(name='protFLEX_prediction',
      version='1.0',
      description='Project that calculates the flexibility of the proteins',
      author='Laura Ciaran Alfano and Neus Pou Amengual',
      author_email='laura.ciaran01@estudiant.upf.edu and neus.pou01@estudiant.upf.edu',
      url='',
      packages = ['protFLEX_prediction'],
      install_requires= REQUIREMENTS,
      scripts=['protFLEX', 'protFLEX_functions', 'protFLEX_graphical_representations', 'dictionaries']
      )

