from setuptools import setup

setup(
    name='aqsoln',
    version='0.6.9',
    author='sizors',
    author_email='sizors@gmail.com',
    packages=[pandas, scipy, sqlite3, molmass],
    url='https://github.com/sizors/aqsoln',
    license='LICENSE.md',
    description='A package for computing solution densities.',
    long_description=open('README.md').read()
)
