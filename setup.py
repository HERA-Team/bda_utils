# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the 2-clause BSD License

"""
Setup file for bda_utils.
"""

import os
import glob
from setuptools import setup, find_packages

setup_args = {
    "name": "bda_utils",
    "author": "The HERA Collaboration",
    "url": "https://github.com/HERA-Team/bda_utils",
    "license": "BSD",
    "description": "a series of tools for supporting BDA",
    "package_dir": {"": "src"},
    "packages": ["bda_utils"],
    "scripts": [fl for fl in glob.glob("scripts/*") if not os.path.isdir(fl)],
    "use_scm_version": True,
    "install_requires": [
        "astropy",
        "numpy",
        "pyuvdata",
        "redis",
    ],
}

if __name__ == "__main__":
    setup(**setup_args)
