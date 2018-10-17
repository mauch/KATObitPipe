#!/usr/bin/env python
from setuptools import setup, find_packages

setup (
    name = "katim",
    version = "trunk",
    description = "",
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
    keywords="kat kat7 ska",
    zip_safe = False,
    # Bitten Test Suite
    #test_suite = "katfile.test.suite",
)
