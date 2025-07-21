# setup.py

from setuptools import setup, find_packages

setup(
    name="splitpea",
    version="0.1.38",
    author="Jeffrey Zhong",
    author_email="jz94@rice.edu",
    description="Splitpea: method for calculating network rewiring changes due to splicing",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ylaboratory/splitpea_package",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'splitpea=splitpea.main:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        "numpy",
        "networkx",
        "intervaltree",
        "ipykernel",
        "matplotlib",
        "plotnine",
        "scikit-learn",
        "adjustText",
        "importlib_resources"
    ],
    include_package_data=True,
    package_data={"splitpea": ["src/*"]},
)
