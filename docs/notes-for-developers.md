# Notes for developers

## Pull Requests

Please use the [Git Flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) model, with a feature branch based off the develop branch, for creating PRs to the [GitHub repo](https://github.com/msk-access/sequence_qc)

## Versioning

To increase the version number use the following command:

`bumpversion (major|minor|patch) --tag`

## Releasing to PyPi and Conda

```text
$ bumpversion (major|minor|patch) --tag
$ python setup.py sdist bdist_wheel
$ twine upload dist/*
$ conda skeleton pypi sequence-qc
---
Optional fixes for potential errors:
    - Resolve ContextualVersionConflict
    - Change "source: url" in meta.yaml to "files.pythonhosted.org/..."
---
$ conda build -c conda-forge -c bioconda sequence-qc
$ anaconda upload /Users/ianjohnson/miniconda3/conda-bld/osx-64/sequence-qc-0.1.12-py37_0.tar.bz2
```

