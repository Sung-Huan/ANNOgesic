from os import path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

setup(
    name='ANNOgesic',
    version='0.7.23',
    packages=['annogesiclib'],
    author='Sung-Huan Yu',
    author_email='sung-huan.yu@uni-wuerzburg.de',
    description='ANNOgesic - A tool for bacterial/archaeal RNA-Seq based genome annotations',
    long_description=long_description,
    url='https://github.com/Sung-Huan/ANNOgesic',
    install_requires=[
        "biopython >= 1.65",
        "matplotlib >= 1.5.0",
        "numpy >= 1.9.2",
        "networkx >= 1.9.1"
    ],
    scripts=['bin/annogesic'],
    license='ISC License (ISCL)',
    classifiers=[
        'License :: OSI Approved :: ISC License (ISCL)',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
