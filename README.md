# SIMBA 2.0

SIMBA 2.0 is an advanced evolutionary bioinformatics tool designed to detect selection pressures in viral genomic data. It provides a significant upgrade from the original SIMBA v1 framework by incorporating more biologically realistic models and modern Bayesian inference.

## Evolution from Version 1

SIMBA 2.0 represents a complete rebuild of the selection detection engine.

### Transition to HKY85

The original version of SIMBA relied on the Nei-Gojobori method. While foundational, the Nei-Gojobori approach assumes that all nucleotide transitions and transversions occur at the same rate. 

SIMBA 2.0 replaces this with the **HKY85 (Hasegawa, Kishino, and Yano 1985)** substitution model. This transition provides several critical improvements:

* **Transition and Transversion Bias** The model explicitly accounts for the fact that certain mutations are more frequent than others.
* **Uneven Base Frequencies** It adjusts for the specific A, C, G, and T distributions found in the SARS-CoV-2 genome.
* **Accuracy** These adjustments result in more precise estimates of synonymous (alpha) and non-synonymous (beta) rates.

### Statistical Improvements

Beyond the substitution model, SIMBA 2.0 introduces several new statistical layers.

* **Empirical Bayesian Grid** Instead of a simple point estimate, the tool uses a Bayesian grid to learn a global prior from the local Zimbabwean dataset. 
* **Numerical Stability via Epsilon Normalisation** This prevents likelihood errors often caused by gaps in lower-quality surveillance sequences.
* **High-Impact Visualisation** The new engine generates publication-standard Manhattan plots with automated domain shading for the RBD and RBM.

## Methodology References

If you use SIMBA 2.0 in your research, please cite the following foundational methods:

1. Hasegawa, M., Kishino, H. and Yano, T. (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA.
2. Murrell, B., et al. (2013). FUBAR: A Fast, Unconstrained Bayesian AppRoximation for Inferring Selection.
3. Nei, M., and Gojobori, T. (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. (Referenced as the baseline for SIMBA v1).

## Usage

Place your genomic sequences in the data directory and run the orchestrator script.

`python run_simba.py`

## Contact

Milton Simbarashe Kambarami | (mskambarami@gmail.com |
 Bioinformatics and Virology Researcher |
 Harare, Zimbabwe
