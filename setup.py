from setuptools import setup, find_packages

setup(
    name="site-annotate",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "Click",
    ],
    entry_points={
        "console_scripts": [
            "site-annotate=site_annotate.cli:main",
        ],
    },
)
