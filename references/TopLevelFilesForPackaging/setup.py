import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PolyMID",
    version="0.0.5",
    author="Nathaniel Vaanti",
    author_email="nv83@cornell.edu",
    description="stable-isotope tracing data analysis software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/VacantiLab/PolyMID",
    packages=setuptools.find_packages(),
    package_data={'': ['*.txt']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

#from: https://stackoverflow.com/questions/1612733/including-non-python-files-with-setup-py
#    package_data is a dict of package names (empty = all packages) to a list of patterns (can include globs). For example, if you want to only specify files within your package, you can do that too:
