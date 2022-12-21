# Supplementary material for "The battle of clean and dirty qubits in the era of partial error correction"
This repository acts as a data repository and codebase to reproduce the results in the paper "The battle of clean and dirty qubits in the era of partial error correction" (https://arxiv.org/abs/2205.13454). 

## Index
### Figure reproduction
The original data to reproduce the numerical figures in the paper is provided in the `results` folder, and code to visualize it can be run in the `Reproduce Figures.ipynb` notebook.
### Data reproduction
The `Recreate data.ipynb` contains the wrapper to recreate the depolarizing noise model of the clean and dirty setup, with code being followable from within the notebook to the relevant `.py` files. A precomputed sample of the output is provided in the `resutls' folder. A huge amount of thanks goes to *Tom O'Leary* for his work on this alternate qiskit code, as the original data was produced with LANL proprietary code. 
