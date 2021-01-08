# Discrete Self Replicating Birth Death Process

Simulates a simple. discrete time stochastic model of a self replicating object, i.e.
```
A -> A + A
A -> 0
```
In addition, objects can mutate, i.e. `A->B`, and B replicates with the same rules.

**Examples**:
- cellular growth and death (mutation rate would e.g. be the emergence of a new genotype)
- PCR amplification

## Installation
**Requirements**: Just a recent python3 installation with pip, e.g. [miniconda](https://docs.conda.io/en/latest/miniconda.html).
```bash
pip install git+https://github.com/redst4r/birthdeathpy.git
```

## Usage
Via the script `birthdeath`:
```bash
birthdeath simulate_uniform_abundance ...
```
This simulates a population of species with equal initial abundance.

```bash
birthdeath simulate_dirichlet_abundance ...
```
Simulates the population of species with initial abundances drawn from a Dirichlet distribution (--alpha parameter!).
This allows for skewed starting distributions.

## Example
```bash
birthdeath simulate_uniform_abundance \
  --n_species 1000 \
  --abundance 1 \
  --rounds 20 \
  --mutation_rate 0.001 \
  --death_probability 0.2 \
  -output_prefix results_
```
This initializes the population with 1000 distinct species, each having one individual (A_i = 1, i=1...1000), replicates for 20 rounds with given death and mutation probabilities (per round) and saves the results under `results_...`.

### Outputs
Outputs are saved under the given `--output_prefix`:
- `...clonesizes.csv`: Sizes of all surviving species after replication
- `...initial_clonesize.png`: Histogram of initial clonesizes
- `...final_clonesize.png`: Histogram of final clonesizes
