# tt+HF analysis macros
This repository contains some useful macros for the tt+HF OSU analysis. To be run at lxplus.

## Get the repo
### Upload with submodules:
```bash
git clone --recursive https://github.com/egorantipov/tt_hf.git
```

### Setup KLFitter
```bash
cd KLFitter
mkdir build && cd build
cmake3 -DBUILTIN_BAT=TRUE ..
make -j
```

## Setup the environment
```bash
source setup.sh
```

## Prepare a histfile
KLFitter implementation is a work in progress. Given the direction on KLFitter setup, one shouldn't run just `roor -l -b prepare_hists_mc.root+`. To prepare a root files with histograms for MC and data run:
```bash
root -l -b load_klf.C+
root -l -b prepare_hists_data.c+
```
Outputs are `hists_mc.root` and `hists_data.root`.

## Draw histograms
To draw histograms prepared by the `prepare_histograms.c` run:
```bash
root -l -b draw_hists.c
```
To draw Data/MC comparison run:
```bash
root -l -b draw_data_mc.c
```
The optut of the macros will be stored in the `Plots/` folder.
