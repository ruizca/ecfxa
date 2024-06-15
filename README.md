# Energy Conversion Factors for X-ray Astronomy

`ecfxa` is a Python package providing tools for the estimation
of ECFs for different astronomical X-ray missions. Currently 
we include XMM-Newton/EPIC, Swift/XRT and eROSITA. ECFs where
calculated consistently for the different missions using the most
updated calibration files at the date of the package release.

X-ray sources detected through photon counter cameras are characterized 
by count-rates at a given energy band. These count-rates (CR) can be 
converted into physical fluxes using the calibrated response of the
detector and assuming an specific spectral model via:

```math
FLUX = CR / ECF
```
where the ECF (Energy Conversion Factor) encapsulate the information
about the response and spectral model for the given energy band. 
For an in-depth explanation of ECF and how they are calculated see
[Mateos et al. 2009](https://doi.org/10.1051/0004-6361/200811409) 
or the [eROSITA ECF tutorial](https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_arfrmf/eROSITA_ECF_tutorial.pdf).

Details about our ECF calculations for the different missions can 
be found in the [corresponding Jupyter notebooks](https://github.com/ruizca/ecfxa/tree/main/calc) included 
in the `ecfxa` GitHub repository.


Installation
------------

    > pip install ecfxa


Examples
--------

    >>> import ecfxa
    >>> ero_ecfs = ecfxa.eROSITA(eband="SOFT")
    >>> ero_ecfs(nh=3e20, gamma=2.0)
    <Quantity 1.18592291e+12 cm2 / erg>

We provide a [Jupyter notebook with examples](https://github.com/ruizca/ecfxa/tree/main/docs/examples.ipynb) 
on how to use `ecfxa` for different X-ray missions. Check the full 
documentation [here](https://ruizca.github.io/ecfxa/ecfxa/ecfs.html).