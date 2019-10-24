# -*- coding: utf-8 -*-
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as fr:
    install_requires = fr.read().splitlines()

setuptools.setup(
    name="geonetworkx",
    version="0.4",
    author="Hugo Chareyre - Artelys",
    author_email="hugo.chareyre@artelys.com",
    description="Python tools for geographic graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ])
