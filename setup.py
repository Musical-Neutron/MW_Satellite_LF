#!/usr/bin/env python3

from setuptools import setup

cmdclass = {}
ext_modules = []

name = 'mw_sats_lf'
version = '1.0'
release = '1.0.0'

setup(
    name=name,
    version=version,
    description='Codes to generate MW satellite galaxy luminosity functions',
    author='Oliver Newton',
    author_email='onewton@cft.edu.pl',
    packages=['mw_sats_lf'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    command_options={},
)
