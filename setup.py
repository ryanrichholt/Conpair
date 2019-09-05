import os
from glob import glob
from setuptools import setup, find_packages


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as fp:
        return fp.read()


def find_package_data():
    files = glob('conpair/data/**', recursive=True)
    # Remove the top-level directory from all the paths
    files = [os.path.sep.join(p.split(os.path.sep)[1:]) for p in files]
    return files


setup(
    name='conpair',
    version='0.2',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    packages=find_packages(exclude=('test',)),
    python_requires='>=3',
    scripts=[
        'scripts/estimate_tumor_normal_contamination.py',
        'scripts/run_gatk_pileup_for_sample.py',
        'scripts/verify_concordance.py'
    ],
    include_package_data=True,
    package_data={
        'conpair': find_package_data()
    },
    install_requires=[
        'scipy'
    ],
)
