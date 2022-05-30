#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "Click>=7.0",
    "numpy",
    "scipy",
    "matplotlib",
    "pandas",
    "seaborn",
    "tqdm",
    "scikit-image",
]

test_requirements = []

setup(
    author="Helium1 LCF Team",
    author_email="helium1bellineq@protonmail.com",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Some tools for the Helium1 team @LCF",
    entry_points={
        "console_scripts": [
            "heliumtools=heliumtools.cli:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="heliumtools",
    name="heliumtools",
    packages=find_packages(include=["heliumtools", "heliumtools.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/quantumatomoptic/heliumtools",
    version="0.1.0",
    zip_safe=False,
)
