# EntropyEstimation

A tool to compute Shannon entropy efficiently. Based on our CAV-22 paper `A Scalable Entropy Estimator` (https://arxiv.org/abs/2206.00921).

## Requirements to run

* Python 2.7+

To install the required libraries, run:

```
python -m pip install -r requirements.txt
```

Additionaly it depends on:

1. [Spur](https://github.com/ZaydH/spur) to sample the satisfying assignment of specification uniformaly at random.
2. [Ganak](https://github.com/meelgroup/ganak) for projected model counting.
3. [WAPS](https://github.com/meelgroup/WAPS): You can choose to sample using WAPS instead of Spur.


In the `dependencies` and `bin` directory, you will find 64-bit x86 Linux compiled binaries for the required dependencies.


## How to Use

```bash
python main.py <input> 
```

This should be able to compute the entropy.

There are different arguments.

### Required
Input CNF file.

### Optional
1. verbose --- different levels of verbosity (0,1,2)
2. epsilon --- a float between 0 and 1. Default 0.8
3. delta --- a float between 0 and 1. Default 0.09
4. sampler --- 1. for Spur and 2. for Waps.
5. counter --- 1. for Ganak
6. timeout --- time out in seconds. Default is 3000s

You can add your custom sampler in `GenerateSamples` and your custom counter in `ComputeCount` subroutine.

### To Run the baseline (enumeration-based approach) to compute the entropy.

```bash
python baseline.py <input>
```

## Benchmarks
Benchmarks are available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6526072.svg)](https://doi.org/10.5281/zenodo.6526072)

## Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new issue](https://github.com/meelgroup/EntropyEstimation/issues). All issues are responded to promptly.

## How to Cite
```
@inproceedings{GJM22,
author={Golia, Priyanka and  Juba, Brendan and  Meel, Kuldeep S.},
title={A Scalable Entropy Estimator},
booktitle={Proceedings of International Conference on Computer-Aided Verification (CAV)},
month={8},
year={2022}
}



```
