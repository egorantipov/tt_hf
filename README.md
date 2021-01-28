# tt+HF analysis macros
This repository contains some useful macros for the tt+HF OSU analysis. To be run at lxplus.

## Get the repo:
```bash
git clone https://github.com/egorantipov/tt_hf.git
```

## Prepare a histfile
To prepare a root file with histograms run:
```bash
prepare_histograms.c+
```
Output of the macro is `hists_mc.root` file.

## Draw histograms
To draw histograms prepared by the `prepare_histograms.c` run:
```bash
draw_histos.c
```
The optut of the macro will be stored in the `Plots/` folder.
