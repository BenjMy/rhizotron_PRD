import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="pyPRD",
    version="0.0.1",
    author="B. Mary",
    author_email="bmary@lbl.gov",
    description="Python package for managing PRD processing, and visualization at the project level",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BenjMy/pyPRD",
    packages=setuptools.find_packages(),
    #packages='pyPRD',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License " + "v3 (LGPLv3)",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=["numpy", "pandas", "pyvista", "resipy",'kneed'], #'icsd_dev',
    dependency_links = ['https://github.com/BenjMy/icsd_dev/tarball/master#egg=icsd_dev-0.1']

) 
