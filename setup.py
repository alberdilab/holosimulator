from setuptools import setup, find_packages

setup(
    name="HoloSimulator",
    version="1.0.0",
    author="Antton Alberdi",
    author_email="antton.alberdi@sund.ku.dk",
    description="Flexible software for simulating holo-omic datasets",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "argparse",
        "PyYAML",
        "requests"
    ],
    entry_points={
        "console_scripts": [
            "holosimulator=holosimulator.cli:main",
        ],
    },
    python_requires=">=3.6",
)
