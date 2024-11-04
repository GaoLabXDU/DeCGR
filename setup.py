
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
    version = "1.1.12",
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
        "pandas",
        "scikit-learn",
        "PyQt5==5.15.11",
        "numpy",
        "wget",
        "cooler",
        "pomegranate==0.15.0",
        "matplotlib",
        "sortedcontainers",
        "scikit-image"
        ]
    )
