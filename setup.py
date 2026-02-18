from setuptools import setup
from osync.__init__ import __version__, __author__, __email__

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="osync",
    version=__version__,
    license="MIT",
    license_file = "LICENSE",
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    description=(
        "osync: an integrated framework for reproducible and efficient processing "
        "of high-throughput sequencing data."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["osync", "osync.workflow.rules", "osync.workflow.utils"],
    package_data={"osync" : ["workflow/*/*"]},
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=[],
    url="https://github.com/carlga/osync",
    keywords=["osync", "multi-omics", "data-driven discovery", "data integration"],
    entry_points={
        "console_scripts": ["osync=osync.main:main"],
    },
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)