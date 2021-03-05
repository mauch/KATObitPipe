#!/usr/bin/env python
from setuptools import setup, find_packages

setup (
    name = "katim",
    version = "0.1.0.dev0",
    description = "Scripts to calibrate MeerKAT data using Obit",
    author = "Tom Mauch",
    author_email = "tmauch@ska.ac.za",
    packages = find_packages(),
    scripts = [
        "scripts/KATContPipe.py",
        "scripts/mvftouvfits.py",
	    "scripts/image_obit.py",
        "scripts/KATCalPipe.py",
	    "scripts/KATZenPipe.py"
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    platforms = [ "OS Independent" ],
    keywords="MeerKAT",
    zip_safe = False,
    install_requires=[
        "astropy",
        "katsdpsigproc",
        "katdal",
        "katpoint",
        "matplotlib",
        "numba",
        "numpy"]
)
