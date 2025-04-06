import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cjson",
    version="0.0.1",
    author="Alexander van Teijlingen",
    author_email="alex@modelmole.com",
    description="A Chemical JSON processing library for Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/beck-vant/cjson",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=["MDAnalysis", "numpy", "ase"],
)
