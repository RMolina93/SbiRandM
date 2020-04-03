#!/usr/bin/env python3
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sbiRandM",
    version="2.0.1",
    py_modules = ['check_stoichometry','exceptions', 'data','superimp','mainSuperimp','mainModeller','pairwise_generator','complex_build','mainWindow','models','pdb_to_fasta'],
    author="Ruben-And-Miguel",
    author_email="ruben.molina-fernandez@upf.edu",
    description="Protein-DNA Complex builder",
    long_description="This is the Structural Bioinformatics project of Ruben Molina and Miguel Dieguez.",
    long_description_content_type="text/markdown",
    url="https://github.com/RMolina93/Structural_Project",
    packages=['sbiRandM'],
    install_requires=['biopython'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
