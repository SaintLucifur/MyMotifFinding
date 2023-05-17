import os
from setuptools import setup, find_packages

curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 1
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'MyMotifFinding/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )
        
setup(
    name='MyMotifFinding',
    version=VERSION,
    description='CSE185 Demo Project',
    author='xxx',
    author_email='xxx',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'mmf=MyMotifFinding.MyMotifFinding:main'
        ],
    },
)