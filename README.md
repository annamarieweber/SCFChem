# SCFChem++
SCFChem++ is a C++ tool for calculating nuclear energy, force, and other properties of molecules using the Self-Consistent Field (SCF) method. This package is designed for use in quantum chemistry and computational materials science, and allows users to easily compute and analyze the electronic structure of molecules.

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

## Usage
To use SCFChem++, run the following command:

```sh
make run mol=molecule
```
Where `molecule` is a file (`molecule.txt`) in the `molecule_files` directory containing the input data for the calculation . The input data should include the atomic coordinates and other necessary information for the calculation. The output will be written to a file in the `output` directory called `molecule.out`. Some example molecules are already included in the molecule_files directory for convenience.

## Docs
[![Publish Docs](https://github.com/annamarieweber/SCFChem/actions/workflows/publish_docs.yml/badge.svg)](https://github.com/annamarieweber/SCFChem/actions/workflows/publish_docs.yml)

## Credits
SCFChem++ was developed and is maintained by Anna Weber.

I hope this readme is helpful! Let me know if you have any questions or need further assistance.
