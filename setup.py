#!/usr/bin/env python3
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sbiRandM",
    version="1.0.0",
    author="Ruben-And-Miguel",
    author_email="ruben.molina-fernandez@upf.edu",
    description="Protein-DNA Complex builder",
    long_description="long_descriptio",
    long_description_content_type="text/markdown",
    url="https://github.com/RMolina93/Structural_Project",
    packages=['sbiRandM','sbiRandM.sbiRandM'],
    install_requires=['biopython'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)