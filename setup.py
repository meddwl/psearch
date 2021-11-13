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
    url="https://github.com/meddwl/psearch",
    packages=setuptools.find_packages(),
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
    install_requires=['pmapper>=0.4.1'],
    entry_points={'console_scripts':
                      ['psearch = psearch.psearch:entry_point',
                       'gen_db = psearch.gen_db:entry_point',
                       'screen_db = psearch.screen_db:entry_point',
                       'prediction = psearch.prediction:entry_point',
                       'external_stat = psearch.scripts.external_statistics:entry_point']},
    include_package_data=True
)
