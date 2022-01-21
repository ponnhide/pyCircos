#! /usr/bin/env python
#
# Copyright (C) Hideto Mori


DESCRIPTION = "circosplot for matplotlib"
LONG_DESCRIPTION = ""

DISTNAME         = 'python-circos'
MAINTAINER       = 'Hideto Mori'
MAINTAINER_EMAIL = 'hidto7592@gmail.com'
URL              = 'https://github.com/ponnhide/pyCircos'
LICENSE          = 'GNU General Public License v3.0'
DOWNLOAD_URL     = 'https://github.com/ponnhide/pyCircos'
VERSION          = '0.2.0'
PYTHON_REQUIRES  = ">=3.7"

INSTALL_REQUIRES = [
    'matplotlib>=3.4',
    'biopython>=1.78',
]

PACKAGES = [
    'pycircos'
]

CLASSIFIERS = [
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
]


if __name__ == "__main__":
    from setuptools import setup
    import sys
    if sys.version_info[:2] < (3, 7):
        raise RuntimeError("pycircos requires python >= 3.7.")

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        packages=PACKAGES,
        classifiers=CLASSIFIERS
    )
