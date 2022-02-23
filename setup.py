from distutils.core import setup
PACKAGES = []
REQUIREMENTS = ['biopython', 'scipy', 'numpy', 'seaborn', 'pandas', 'matplotlib']

setup(name='Laura_Neus',
      version='1.0',
      description='Project that calculates the flexibility of the proteins',
      author='Laura Ciaran Alfano and Neus Pou Amengual',
      author_email='laura.ciaran01@estudiant.upf.edu and neus.pou01@estudiant.upf.edu',
      url='?',
      packages = PACKAGES,
      requires= REQUIREMENTS,
     )