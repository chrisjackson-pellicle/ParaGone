#!/usr/bin/env python

import setuptools

resolve_paralogs_scripts = ['paragone']
resolve_paralogs_description = 'Paralogy resolution using Yang and Smith algorithms'
resolve_paralogs_url = 'https://github.com/chrisjackson-pellicle/yang_and_smith_paralogy_resolution'
resolve_paralogs_entry_points = {'console_scripts': ['paragone = paragone.paragone_main:main']}

setuptools.setup(name='paragone',
                 version='0.0.5rc',
                 packages=setuptools.find_packages(),
                 author='Chris Jackson',
                 author_email='chris.jackson@rbg.vic.gov.au',
                 description=resolve_paralogs_description,
                 keywords='paralog resolution, orthology inference',
                 url=resolve_paralogs_url,
                 entry_points=resolve_paralogs_entry_points)
