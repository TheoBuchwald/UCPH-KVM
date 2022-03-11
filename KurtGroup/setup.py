from setuptools import setup


setup(
    name='KurtGroup',
    version='1.1.11a0',
    author='Theo Juncker von Buchwald',
    author_email='fnc970@alumni.ku.dk',
    packages=['Kurt'],
    package_data={'Kurt': ['*.py', '*.json']},
    license='LICENCSE',
    description='A package containing scripts and programs developed and used in Professor Kurt V. Mikkelsens group at the University of Copenhagen.',
    long_description=open('README.md').read(),
    install_requires=['numpy >= 1.13.0', 'requests >= 2.4', 'ase >= 3.19.0', 'matplotlib'],
    python_requires='>=3.6'
)