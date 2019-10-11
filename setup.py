# coding=utf-8

import setuptools

def read_version_string(filename):
    versionstring = None
    with open(filename) as h:
        for line in h:
            try:
                identifier, equalsign, versionstring = line.split()
                if identifier == '__version__' and equalsign == '=':
                    return versionstring.strip("'")
            except:
                continue        # Just read past lines to fitting the pattern
    if versionstring is None:
        return 'unknown_version'

__version__ = read_version_string('modelestimator/version.py')


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='modelestimator-v2',
    version=__version__,
    author='Ruben Ridderström and Lars Arvestad',
    author_email='ruben.ridderstrom@gmail.com',
    description='Program for estimating amino acid replacement rates',
    long_description=long_description,
    long_description_content_type="text/markdown; charset=UTF-8",
    url='https://github.com/arvestad/modelestimator-v2',
    license='GPLv3',
    packages = setuptools.find_packages(),
    entry_points={
        'console_scripts':[
            'modelestimator = modelestimator.main:main'
        ]
    },
    install_requires=[
        'argparse',
        "scipy",
        "numpy",
        "biopython"
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest']
)
