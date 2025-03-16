# setup.py

from setuptools import setup, find_packages

setup(
    name="splitpea",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A short description of what splitpea does",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/raissinging/splitpea_package",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'splitpea=splitpea.main:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        "numpy==1.24.3",
        "networkx==3.1",
        "intervaltree==3.1.0",
        "ipykernel",
        "karateclub",
        "pandas",
        "matplotlib",
        "plotnine",
        "scikit-learn",
        "mygene"
    ],
)
