#  Online supplement to "A genetic algorithm for fair land allocation"

This repository contains instances, experimental data, and code for the paper Gliesch, Ritt, Moreira, [A genetic algorithm for fair land allocation, Proc. GECCO 2017, 793-800](http://dx.doi.org/10.1145/3071178.3071313).

## Tables

In subdirectory `data` you can find the following supplementary data to the results presented in the paper:

  * Table 1: [Information on the instances](data/instance-stats.csv).
  * Table 2: [irace parameter files and results](irace). 
  * Table 3: results of the calibration of the optimal batch size: [fixed time](data/batch-sizes.csv) and [50 replications](data/batch-sizes-50-generations.csv). 
  * Tables 4 and 5: scalability and effectiveness experiments: [constructive algorithm](data/results-constructive-30-min.csv), [BFS](data/results-naive-30-min.csv) and [GA](data/results-genetic-30-min.csv). 
  * Table 6: comparison of [our GA](data/results-genetic-30-min.csv) and the [manual allocation](data/results-incra.csv). 

## Instance generator

The instance generator is in the subdirectory `instances/instance-generator`. You need to Boost libraries to compile. If they're installed, use
```bash
cd instances/instance-generator; make -j
```
to build it.

## Source code

The code is contained in the subdirectory `src`. You need to Boost libraries to compile. To build:
```bash
cd src; make -j
```	

## How to cite

```bibtex
@InProceedings{Gliesch.etal/2017,
  author =    {Alex Gliesch and Marcus Ritt and Mayron C. O. Moreira},
  title =     {A genetic algorithm for fair land allocation},
  booktitle = {Proc. 19th Conf. Genetic Evol. Comput.},
  editor =    {Peter A. Bosnan},
  crossref =  {gecco2017},
  year =      {2017},
  pages =     {793--800},
  address =   {Berlin},
  publisher = {ACM Press},
  doi =       {10.1145/3071178.3071313},
  isbn =      {9781450349208}
}
```
