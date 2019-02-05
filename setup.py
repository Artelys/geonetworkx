"""
    File name: setup.py
    Author: Artelys
    Creation date: 15/01/2019
    Python Version: 3.6
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = ["pyproj>=1.9.6",
                    "geopy>=1.12.0",
                    "geopandas>=0.4.0",
                    "networkx>=2.2",
                    "numpy>=1.15.4",
                    "shapely>=1.2.18",
                    "scipy>=1.0.1"]

setuptools.setup(
    name="geonetworkx",
    version="0.1",
    author="Artelys",
    author_email="hugo.chareyre@artelys.com",
    description="Python tools for geographic graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ])