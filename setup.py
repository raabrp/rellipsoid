import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rellipsoid",
    version="0.0.1",
    author="Reilly Raab",
    author_email="raabrp@gmail.com",
    description="Ellipsoidal model of planet with WSG84 values for Earth",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/raabrp/rellipsoid",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
