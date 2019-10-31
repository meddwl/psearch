import setuptools
from os import path
import psearch

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="psearch",
    version=psearch.__version__,
    author="Alina Kutlushina, Pavel Polishchuk",
    author_email="alina.kutlushina@pharminnotech.com, pavel_polishchuk@ukr.net",
    description="PSearch: ligand-based pharmacophore modeling and screening",
    long_description=long_description,
    long_description_content_type="text/markdown",
#######    url="https://github.com/DrrDom/pmapper",
    packages=['psearch'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
    extras_require={
    },
    install_require={
        'pmapper': ['pmapper>=0.2'],
    },
    entry_points={'console_scripts':
                      ['get_pharm_counts = pmapper.get_pharm_counts:entry_point']},
    include_package_data=True
)
