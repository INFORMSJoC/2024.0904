[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An LBBD for the CRSP with A PBS

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[A Logic-Based Benders Decomposition for the Car Resequencing Problem with a Painted Body Storage](https://doi.org/10.1287/ijoc.2024.0904) by X. Guo, J.-F. Côté, C. Zhang, and L. Miao.


## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0904

https://doi.org/10.1287/ijoc.2024.0904.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{guo2025lbbd,
  author =        {Guo, Xinyi and C{\^o}t{\'e}, Jean-Fran{\c{c}}ois and Zhang, Canrong and Miao, Lixin},
  publisher =     {INFORMS Journal on Computing},
  title =         {{A Logic-Based Benders Decomposition for the Car Resequencing Problem with a Painted Body Storage}},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0904.cd},
  url =           {https://github.com/INFORMSJoC/2024.0904},
  note =          {Available for download at https://github.com/INFORMSJoC/2024.0904},
}
```

## Description

The purpose of this repository is to enable the reproduction of the tables and figures presented in the paper and its Supplementary Material.


## Instances

All instances used in the paper are included in folder `data`, in which

- `BasicData`: includes all instances for Algorithm Comparison
  
  * `50_1.txt` in folder `RealCase` is the 1*st* real-world instance with 50 cars
  
- `ExtendedData`: includes all instances for Sensitivity Analysis

  * all instances in folder `233-B=2 M=3 J=3` use 2 _bodies_, 3 _colors_, and 3 _configurations_

 
## Algorithms

The codes of all algorithms in the paper are available in folder `src`, in which

- `Complete MIP` includes codes for solving the complete formulation by Gurobi.
- `kSSP-G1` includes codes for executing the NLBBD based on G1.
- `NLBBD & NLBBD_G2` includes codes for executing the basic NLBBD and the one based on G2.
- `Standard LBBD` includes codes for executing the standard two-level LBBD.

## Results

Detailed results for the tables and figures in the paper are included in the folder `results`.


## Computing Environment

The C++ codes and Bash scripts are intended for use in a Linux environment equipped with standard libraries and utilities such as the g++ compiler and the Bash shell.
