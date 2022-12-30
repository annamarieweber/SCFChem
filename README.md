# SCFChem
SCFChem is a C++ tool for calculating nuclear energy, force, and other properties of molecules using the Self-Consistent Field (SCF) method. This package is designed for use in quantum chemistry and computational materials science, and allows users to easily compute and analyze the electronic structure of molecules.

## Installation

### To install SCFChem++, follow these steps:

Clone the repository using the following command:
```sh
git clone https://github.com/your-username/SCFChem++.git
```
Navigate to the directory where you cloned the repository:

```sh
cd SCFChem++
```

Compile the code using the provided makefile:

```sh
make build all
```

### Dependencies
SCFChem requires the armadillo library to be installed on your system. Armadillo can be downloaded from https://arma.sourceforge.net/download.html. You should be able to find installation instructions there as well. 

If your armadillo header and library files are installed in a location other than `/usr/local/include` and `/usr/local/lib` respectively be sure to set the environment variable `LOCAL_PATHS` using one of the methods described below.

#### Set LOCAL_PATHS at the directory level with direnv
If you have direnv set up on your system you can set the environment variable by creating a .envrc file and adding the following line to it: ```LOCAL_PATHS="-I/path/to/armadillo/headers -L/path/to/armadillo/libraries"``` and then running direnv allow. Alternatively the following line will also complete these steps for you `echo LOCAL_PATHS="-I/path/to/armadillo/headers -L/path/to/armadillo/libraries"`.


#### Set LOCAL_PATHS in your shell profile
Add ```export LOCAL_PATHS="-I/path/to/armadillo/headers -L/path/to/armadillo/libraries"``` to the profile file for your shell of choice. For example, if you use zsh add the line to your `~/.zshrc` file.

#### Temporarily set LOCAL_PATHS in terminal
Run ```export LOCAL_PATHS="-I/path/to/armadillo/headers -L/path/to/armadillo/libraries"``` in the terminal prior to running make commands. note that the environment variable will only be set for that session and you will need to re run the command each time you start a new session. 

`path/to/armadillo/headers` and `/path/to/armadillo/libraries` should be the parent directories where these items are stored. For example if your armadillo installation has the header files at `/usr/share/include/armadillo` and `/usr/share/include/armadillo_bits` and library files at `/usr/share/lib/libarmadillo.dylib` then your `LOCAL_PATHS` environment variable should be `"-I/usr/share/include -L/usr/share/lib"`

## Usage
To use SCFChem, run the following command:

```sh
make run mol=molecule
```
Where `molecule` is a file (`molecule.txt`) in the `molecule_files` directory containing the input data for the calculation . The input data should include the atomic coordinates and other necessary information for the calculation. The output will be written to a file in the `output` directory called `molecule.out`. Some example molecules are already included in the molecule_files directory for convenience.

## Docs
[![Publish Docs](https://github.com/annamarieweber/SCFChem/actions/workflows/publish_docs.yml/badge.svg)](https://github.com/annamarieweber/SCFChem/actions/workflows/publish_docs.yml)

## Credits
SCFChem was developed and is maintained by Anna Weber.

I hope this readme is helpful! Let me know if you have any questions or need further assistance.
