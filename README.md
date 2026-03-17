# SIMBA 2.0

SIMBA 2.0 is an evolutionary bioinformatics pipeline designed to detect selection pressures in viral genomic data. It specifically focuses on identifying codon sites under positive selection within the SARS-CoV-2 Spike protein. 

The tool bridges the gap between complex statistical inference and biological interpretation by automatically mapping evolutionary hotspots to functional protein domains.

## Core Methodology

The pipeline employs an empirical Bayesian framework to estimate synonymous and non-synonymous substitution rates.

* **HKY85 Engine** The model accounts for transition and transversion bias alongside uneven nucleotide frequencies.
* **Bayesian Grid Inference** Site-specific rates are estimated using a global prior learned from the entire dataset.
* **Numerical Stability** Integrated epsilon-normalisation prevents calculation errors in regions with low sequence coverage.

## Installation

Ensure you have a Python environment ready. You may install the necessary dependencies using the provided requirements file.

`pip install -r requirements.txt`

## Getting Started

1. Place your genomic sequences in the data directory.
2. Ensure your reference sequence is correctly configured in the orchestrator script.
3. Run the primary script from your terminal.

`python run_simba.py`

## Visualisation

SIMBA 2.0 generates high-resolution graphics suitable for peer-reviewed publications. 

* **Selection Landscape** A Manhattan-style plot showing the probability of selection across the Spike gene with domain-specific shading.
* **Rate Divergence** A scatter plot illustrating the relationship between alpha and beta rates to highlight evolutionary hotspots.

## Contact and Contributions

If you are interested in collaborating on viral evolution research or localising these tools for different genomic contexts, feel free to reach out.