# Helicity Computation

This script is intended for the computation of the helicity of crystalline structures. A full description of the method, including its strengths and weaknesses, and its application to real systems can be found in [Phys. Rev. B 110, 174112 (2024)].
The work was developped by Fernando Gómez-Ortiz at the University of Liège working under the supervision of Eric Bousquet, PI of the CHRYSALID PDR project No. 40003544.

The work was developped at the University of Liège and was funded by the Fonds de la Recherche Scientifique (FNRS) through the PDR project CHRYSALID No.40003544 and the Consortium des Équipements de Calcul Intensif (CÉCI), funded by the F.R.S.-FNRS under Grant No. 2.5020.11
## Input

The script takes three parameters:

1. **a (int)**: The number of atoms in the home unit cell of the materials.
2. **b (string)**: The path to a file containing the atomic positions in fractional units of the reference high symmetry structure in the space group of the low symmetry phase.
3. **c (string)**: The path to a file containing the atomic positions in fractional units of the low symmetry phase, in the same order as the previous file (one-to-one mapping).

## Output

The script returns the value of the helicity for the crystalline phase:

- **float**: The helicity of the structure.

## Install

pip install Helicity_Computation_Crystals

## Usage

To run the script, use the following command:

```bash
helicity_computation 24 examples/ref_hs examples/chir_aa
```

Or

```bash
helicity_computation 24 examples/ref_hs examples/ls_a0
```
