project: fortperc
project_dir: ./src
output_dir: ./docs
author: Anders Johansson
project_github: https://github.com/anjohan/fortperc
summary: fortperc - A percolation library in modern Fortran
github: https://github.com/anjohan
graph: true

## Introduction
This library consists of two main modules:

* `hk`, which contains an implementation of the Hoshen-Kopelman algorithm.
* `percolation`, which contains procedures for doing various percolation related measurements, such as the probability for having a percolating cluster (one which allows transport across the system).

Additionally, there is `utilities`, a small module containing utility procedures such as `linspace` and `linfit`.
