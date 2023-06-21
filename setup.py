# copied from Dr. Melissa Gymrek with minor modifications
# https://github.com/gymreklab/cse185-demo-project/blob/main/setup.py

import os
from setuptools import setup, find_packages

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 2
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'inch/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )


setup(
    name='inch',
    version=VERSION,
    description='CSE185 Class Project',
    author='Faith Okamoto',
    author_email='fokamoto@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "inch=inch.inch:main"
        ],
    },
)