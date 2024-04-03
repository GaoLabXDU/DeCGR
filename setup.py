
"""
Setup script for DeCGR.

"""
import os, sys, glob
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major!=3) or (sys.version_info.minor<7):
    print('PYTHON 3.7+ IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)

setuptools.setup(
    name = 'DeCGR',
    version = "1.0.5",
    author = "Li Junping",
    author_email = 'lijunping02@qq.com',
    url = 'https://github.com/GaoLabXDU/DeCGR',
    description = 'An interactive toolkit for deciphering CGRs from Hi-C data',
    keywords = ("Hi-C", "complex rearrangements", "3D genome"),
    scripts = glob.glob('scripts/*'),
    packages = setuptools.find_packages(),
    include_package_data = True,
    platforms = "any",
    license="MIT Licence",
    install_requires = [
        "pandas==2.2.0",
        "pyBigWig==0.3.22",
        "scikit-learn== 1.4.0",
        "PyQt5==5.15.10",
        "numpy==1.26.3",
        "iced==0.5.13",
        "cooler==0.9.3",
        "pomegranate==0.14.8",
        "matplotlib==3.8.2",
        "sortedcontainers==2.4.0",
        ]
    )
