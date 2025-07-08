from setuptools import setup, find_packages

setup(
    name="amplocator",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "tensorflow",
        "biopython",
        "pandas",
        "numpy"
    ],
    entry_points={
        "console_scripts": [
            "amplocator=amplocator.cli:main"
        ]
    }
)
