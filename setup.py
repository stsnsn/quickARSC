from setuptools import setup, find_packages

setup(
    name="arsc",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["biopython>=1.79"],
    entry_points={
        "console_scripts": [
            "ARSC=ARSC.main:main",
        ],
    },
    author="Satoshi Nishino",
    author_email="satoshi-nishino@g.ecc.u-tokyo.ac.jp",
    description="Compute ARSC (N/C/S/AvgResMW) from protein fasta files",
    url="https://github.com/stsnsn/ARSC",
    license="GPL-2.0",
)