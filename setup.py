
from setuptools import setup

setup(
    name='Kurt',
    version='0.0.1',
    author='Theo Juncker von Buchwald',
    author_email='fnc970@alumni.ku.dk',
    packages=['dependencies'],
    license='LICENCSE',
    description='A package containing scripts and programs developed and used in Professor Kurt V. Mikkelsens group at the University of Copenhagen.',
    long_description=open('README.md').read(),
    install_requires=['numpy >= 1.18.0'],
    python_requires='>=3.6'
)