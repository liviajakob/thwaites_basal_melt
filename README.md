# Thwaites basal melt
Public repository for generating basal melt rates over Thwaites.

## Installation
To start with, please ensure that you have some form of [Anaconda](https://www.anaconda.com/products/distribution)
installed (for Python 3.9). Earthwave prefers to use 
[Miniconda](https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh) on Linux.

After cloning the repository, change directories into the newly created folder.

Next, create a new conda environment:

```
conda create --name my_env python=3.9.7
conda activate my_env
```

Note that Earthwave uses VSCode (v1.69.1 at time of writing) as our integrated development environment. Earthwave's standard VSCode environment settings are included with this repository. If you wish to use VSCode, you'll need to install the Python and Pylance plugins, and [tell VSCode to use the interpreter within the environment you've just created](https://code.visualstudio.com/docs/python/environments#_select-and-activate-an-environment).

For a development install, ensure that you install all of the automated testing and linting requirements:
```
pip install -r .github/test_requirements.txt
```

Then install this package for development:
```
pip install -e .
```

You should now be able to run and pass the unit tests simply by running:
```
pytest
```

