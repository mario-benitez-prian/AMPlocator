from setuptools import setup, find_packages

setup(
    name="amplocator",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "tensorflow-cpu>=2.11.0",  # Estable y compatible con la mayoría de CPUs modernas
        "biopython>=1.79",         # Compatible con Python 3.8+ y muy estable
        "pandas>=1.3.0",           # Buen soporte de funciones, sin romper compatibilidad
        "numpy>=1.20.0"            # Compatible con TensorFlow 2.11+
    ],
    entry_points={
        "console_scripts": [
            "amplocator=amplocator.cli:main"
        ]
    }
)
