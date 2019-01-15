"""
    File name: setup.py
    Author: Artelys
    Creation date: 15/01/2019
    Python Version: 3.6
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="geonetworkx",
    version="0.1",
    author="Artelys",
    author_email="hugo.chareyre@artelys.com",
    description="Python tools for geographic graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ])