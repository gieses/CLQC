from distutils.core import setup

setup(
    name='CLQC',
    version='1.1',
    author='Sven H. Giese',
    author_email='sven.giese@tu-berlin.de',
    packages=['CLQC'],
    scripts=['bin/CLQC_ContactMap.py','bin/CLQC_Distogram.py', 'bin/CLQC_CoverageProfile.py'],
    url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE.txt',
    description='Package for fast quality control for CLMS data.',
    long_description=open('README.txt').read(),
    install_requires =["HTSeq >= 0.0.0"],
)
