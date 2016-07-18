import os
from os.path import relpath, join
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files

setup(
    name = "openforcefield",
    version = "0.1.0",
    author = "John D. Chodera",
    author_email = "john.chodera@choderalab.org",
    description = ("Open tools for Bayesian forcefield parameterization"),
    license = "GNU Lesser General Public License (LGPL), Version 3",
    keywords = "Bayesian inference; forcefield parmaeterization; ThermoML",
    url = "http://github.com/open-forcefield-group/open-forcefield-tools",
    packages=['openforcefield', 'openforcefield/tests', 'openforcefield/data'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU Lesser General Public License (LGPL), Version 3",
    ],
    #entry_points={'console_scripts': ['smarty = smarty.cli:main']},
    package_data={'openforcefield': find_package_data('openforcefield/data', 'openforcefield')},
)
