# Cooling of the CGM under shock, adiabatic and continuous compression

This repository contains the code for *[paper title]* (Wu, in prep.).

## Layout

```
reproduce_figures.ipynb   Figures 1-5 of the paper, in order
cooling.py                constants, Wiersma et al. (2009) cooling table, tcool_P
shock_transition.py       shock_tran, time_diff                   (Fig. 2)
adiabatic_process.py      shock_tran, adi_tran, time_diff, kai_func (Figs. 1, 3)
infall_galaxy.py          cont_radi, cool_radiation, entropy_diff,
                          cool_radiation_single                    (Figs. 4, 5)
data/                     cooling table + Figure 5 grids
exploration/              original working notebooks
```

## Usage

```bash
pip install -r requirements.txt
# put z_0.000.hdf5 (Wiersma et al. 2009, z = 0) in data/
jupyter notebook reproduce_figures.ipynb
```

The Figure 5 grids (3 × 200 × 200 ODE solves) take a while to compute the first time. After that they are cached in `data/`.

## Citation

Please cite the paper (reference to be added) and Wiersma, Schaye & Smith (2009), MNRAS 393, 99.
