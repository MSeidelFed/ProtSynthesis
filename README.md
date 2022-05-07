# ProtSynthesis
R package to calculate fractional protein synthesis rates

For detailed usage please check out the [isotopeEnrichment](https://github.com/mgleeming/isotopeEnrichment/blob/master/README.md) repository.

## Installation

```
library(devtools)

### Get the latest installation of RandoDiStats (KineticMSI depends on it)

devtools::install_github("MSeidelFed/RandodiStats_package")
library(RandoDiStats)

### Install KineticMSI

devtools::install_github("MSeidelFed/ProtSynthesis")
library(ProtSynthesis)

```

## Getting the exemplary datasets from the system directory after package installation

These datasets were used to produce the results in XX journal publication.

```
system.file("extdata", package = "ProtSynthesis")
```
