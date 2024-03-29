#!/usr/bin/env python
import os
from setuptools import setup, find_packages


with open('docs/README.md') as readme_file:
    readme = readme_file.read()


def req_file(filename):
    """
    We're using a requirements.txt file so that pyup.io can use this for security checks

    :param filename:
    :return str:
    """
    with open(filename) as f:
        content = f.readlines()
        content = filter(lambda x: not x.startswith("#"), content)
    return [x.strip() for x in content]


def get_package_files(file_type):
    """
    helper function to recursively extract specific file types from the repository.

    :param directory, file_type:
    :return str:
    """
    paths = []
    for (path, directories, filenames) in os.walk(
        os.path.dirname(os.path.abspath(__file__))
    ):
        for filename in filenames:
            if not filename.endswith(file_type):
                continue
            paths.append(os.path.join("..", path, filename))
    return paths


with open(os.path.join(os.path.dirname(__file__), "sequence_qc/VERSION"), "r") as fh:
    __version__ = fh.read().strip()


setup(
    author="Ian Johnson",
    author_email='ionox0@gmail.com',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Package for doing various ad-hoc quality control steps from MSK-ACCESS generated FASTQ or BAM files",
    entry_points={
        'console_scripts': [
            'calculate_noise=sequence_qc.cli:calculate_noise',
        ],
    },
    install_requires=req_file("requirements.txt"),
    license="Apache Software License 2.0",
    long_description=readme,
    include_package_data=True,
    keywords='sequence_qc',
    name='sequence_qc',
    packages=find_packages(include=['sequence_qc', 'sequence_qc.*']),
    package_data={
        "": ['requirements.txt', 'requirements_dev.txt'],
    },
    test_suite='tests',
    url='https://github.com/msk-access/sequence_qc',
    version=__version__,
    zip_safe=False,
)
