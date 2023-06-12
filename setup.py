from setuptools import setup, find_packages

setup(
    name = "photozpy",
    version = "0.1",
    author = "Yong Sheng",
    author_email = "sheng2@clemson.edu",
    description = "Automatic pipeline for data analysis UVOT and SARA images.",

    packages = find_packages('src', exclude = ["tests"]),

    install_requires(numpy,
                     pandas,
                     tqdm,
                     astropy,
                     swifttools)

)
