# ProtSynthesis
R package to calculate fractional protein synthesis rates

For detailed usage of this package and to understand how to produce the input data for its functions please check out the [isotopeEnrichment](https://github.com/mgleeming/isotopeEnrichment/blob/master/README.md) repository, which contains the python dependencies on which this R package relies.

For documentation on the functions provided and their parameters, please consult the provided [PDF file](https://github.com/MSeidelFed/ProtSynthesis/blob/main/ProtSynthesis_0.1.0.pdf).

## Installation

```
library(devtools)

### Get the latest installation of RandoDiStats (ProtSynthesis depends on it)

devtools::install_github("MSeidelFed/RandodiStats_package")
library(RandoDiStats)

### Install ProtSynthesis

devtools::install_github("MSeidelFed/ProtSynthesis")
library(ProtSynthesis)

```

## Getting the exemplary datasets from the system directory after package installation

These datasets were used to produce the results in XX journal publication.

```
system.file("extdata", package = "ProtSynthesis")
```
