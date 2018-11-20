import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FrAG-PELE",
    version="1.0.0",
    author="Carles Perez Lopez & Dani Soler",
    author_email="carles.perez@bsc.es",
    description="FrAG is a bioinformatic tool focused on fragment-based drug design.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/danielSoler93/Ligand_growing.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.4",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
)
