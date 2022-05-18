import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AJMAnalysis", 
    version="0.9",
    author="Rastko Sknepnek",
    author_email="r.sknepnek@dundee.ac.uk",
    description="Analysis toolking for the AJM package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sknepneklab/AJM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: to be determined",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.4',
)
