# coding=utf-8

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('modelestimator/version.py') as fh:
    exec(fh.read())

setuptools.setup(
	name='modelestimator-v2',
	version='0.0.6',
	author='Ruben Ridderström and Lars Arvestad',
	author_email='ruben.ridderstrom@gmail.com',
	description='Program for estimating amino acid replacement rates',
	long_description=long_description,
	long_description_content_type="text/plain; charset=UTF-8",
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
