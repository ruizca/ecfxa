# -*- coding: utf-8 -*-
import gzip
import json
from datetime import datetime
from importlib import resources

import numpy as np
from astropy import units as u
from astropy.time import Time
from scipy.interpolate import RectBivariateSpline

from .singleton import SingletonMeta


class eROSITA:
    """
    Energy Conversion Factors (ECFs) for the eROSITA instrument on-board the Spektr-RG mission.

    We include ECF estimates for the different energy bands used in eROSITA products. ECFs 
    were estimated assuming an absorbed powerlaw. Once an eROSITA instance is initialized for 
    a given energy band, ECFs values can be obtained for different values of NH and photon 
    index. Returned ECFs are in units of counts × cm² / erg and it is possible to include 
    the correction due to absorption.
    
    Examples
    --------
    - eROSITA ECF at energy band P3 (1-2 keV), assuming an absorbed powerlaw with a 
    Hydrogen column density of 5×10²¹ cm-2 and photon index 1.9:

    >>> ero_ecf = eROSITA(eband="P3")

    ECF with no absorption correction:
    
    >>> ero_ecf(nh=5e21, gamma=1.9)
    <Quantity 9.70689828e+11 cm2 / erg>
    
    ECF including absorption correction:

    >>> ero_ecf(5e21, 1.9, abscorr=True)
    <Quantity 5.67143189e+11 cm2 / erg>
     
    - Show available energy bands:
    
    >>> eROSITA.ebands
    """

    ebands = {
        ## Bands from ERASS1 catalogue
        # Broad
        "1": (0.2, 2.3),
        "5": (0.5, 2.0),
        "P1": (0.2, 0.5),
        "P2": (0.5, 1.0),
        "P3": (1.0, 2.0),
        "P4": (2.0, 5.0),
        "P5": (5.0, 8.0),
        "P6": (4.0, 10.0),
        # Narrow
        "P7": (5.1, 6.1),
        "P8": (6.2, 7.1),
        "P9": (7.2, 8.2),
        ## From upper-limit server (Tubin-Arenas+2024)
        "021": (0.2, 0.6),
        "022": (0.6, 2.3),
        "023": (2.3, 5.0),
        "02e": (0.2, 5.0),
        ## Standard
        "SOFT": (0.5, 2.0),
        "HARD": (2.0, 10.0),
    }

    epochs = {
        "e1": ("2019-10-17", Time(datetime.now())),
    }

    def __init__(self, eband="SOFT", date=None):
        self._ecf = eROSITAECFValues()

        # Attributes set in the parse method
        self._parse_args(eband, date)
        self._interpolators = self._set_interpolators()

    def _parse_args(self, eband, date):
        self.eband = self._parse_eband(eband)
        self.epoch = self._parse_date(date)
        
    def _parse_eband(self, eband: str) -> str:
        if eband not in self.ebands:
            raise ValueError(f"Unknown energy band: {eband}")
        
        return eband
        
    def _parse_date(self, date) -> str:
        if date is None:
            epoch = "e1"

        else:
            epoch = None
            date = Time(date)

            for key_epoch, values in self.epochs.items():
                date_min = Time(values[0])
                date_max = Time(values[1])

                if date_min <= date <= date_max:
                    epoch = key_epoch
                    break

        if epoch is None:
            raise ValueError("Date is not compatible with the eROSITA mission.")
        
        return epoch

    def _set_interpolators(self):
        ecf_values_nocorr = np.array(
            self._ecf.nocorr[self.epoch][self.eband]
        )
        ecf_values_abscorr = np.array(
            self._ecf.abscorr[self.epoch][self.eband]
        )
        
        interpolator = {
            "nocorr": RectBivariateSpline(
                self._ecf.nocorr["lognh"],
                self._ecf.nocorr["gamma"],
                ecf_values_nocorr,
                kx=1,
                ky=1,
            ),
            "abscorr": RectBivariateSpline(
                self._ecf.abscorr["lognh"],
                self._ecf.abscorr["gamma"],
                ecf_values_abscorr,
                kx=1,
                ky=1,
            ), 
        }

        return interpolator

    def __call__(self, nh=3e20, gamma=1.7, abscorr=False):
        lognh = np.log10(nh)

        # Keep values of lognh and gamma between interpolation limits
        lognh = np.maximum(lognh, self._ecf.nocorr["lognh"][0])
        lognh = np.minimum(lognh, self._ecf.nocorr["lognh"][-1])

        gamma = np.maximum(gamma, self._ecf.nocorr["gamma"][0])
        gamma = np.minimum(gamma, self._ecf.nocorr["gamma"][-1])

        if abscorr:
            ecf = self._interpolators["abscorr"].ev(lognh, gamma)
        else:
            ecf = self._interpolators["nocorr"].ev(lognh, gamma)

        return ecf * 1e11 << u.erg**-1 * u.cm**2


class eROSITAECFValues(metaclass=SingletonMeta):
    data_path = resources.files("ecfxa.data")

    def __init__(self) -> None:
        # Values with no Galactic absorption correction
        with resources.as_file(self.data_path / "erosita_ecfs.json.gz") as ecf_file:
            with gzip.open(ecf_file) as fp:
                self.nocorr = json.load(fp)

        # Values taking into account Galactic absorption correction
        with resources.as_file(self.data_path / "erosita_abscorr_ecfs.json.gz") as ecf_file:
            with gzip.open(ecf_file) as fp:
                self.abscorr = json.load(fp)
