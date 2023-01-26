# Supplementary material for "The battle of clean and dirty qubits in the era of partial error correction"
This repository acts as a data repository and codebase to reproduce the results in the paper "The battle of clean and dirty qubits in the era of partial error correction" (https://arxiv.org/abs/2205.13454). 

## Index
### Figure reproduction
The original data to reproduce the numerical figures in the paper is provided in the `results` folder, and code to visualize it can be run in the `Reproduce Figures.ipynb` notebook. The saved figures can be found in the `plots` folder. 
### Data reproduction
The `Recreate data.ipynb` contains the wrapper to recreate the depolarizing noise model of the clean and dirty setup, with code being followable from within the notebook to the relevant `.py` files. A precomputed sample of the output is provided in the `results` folder. I would like to thank Tom O'Leary for his work on the qiskit implementation. 

### Requirements
Pandas (1.3.3), seaborn (0.11.2), numpy (1.23.4), qiskit (0.22.1), mthree (1.1.0.dev0), matplotlib (3.4.3) and pyarrow (8.0.0).
Alternate versions will work, but these are guaranteed to function. 


If you find this work useful, please consider citing the original paper:
```
@misc{https://doi.org/10.48550/arxiv.2205.13454,
  doi = {10.48550/ARXIV.2205.13454},
  url = {https://arxiv.org/abs/2205.13454},
  author = {Bultrini, Daniel and Wang, Samson and Czarnik, Piotr and Gordon, Max Hunter and Cerezo, M. and Coles, Patrick J. and Cincio, Lukasz},
  keywords = {Quantum Physics (quant-ph), FOS: Physical sciences, FOS: Physical sciences},
  title = {The battle of clean and dirty qubits in the era of partial error correction},
  publisher = {arXiv},
  year = {2022},
  copyright = {arXiv.org perpetual, non-exclusive license}
}
```
