# -*- coding: utf-8 -*-
import gzip
import json
from datetime import datetime
from enum import Enum
from importlib import resources

import numpy as np
from astropy import units as u
from astropy.time import Time
from scipy.interpolate import RectBivariateSpline

from .singleton import SingletonMeta


class XMMEPIC:
    """
    Energy Conversion Factors (ECFs) for the EPIC cameras on-board the XMM-Newton observatory.

    We include ECF estimates for the three detectors of XMM-EPIC (EPN, EMOS1 and EMOS2),
    their corresponding filters and operation modes, and the different epochs defined by the 
    XMM-Newton calibration team. ECFs were estimated assuming an absorbed powerlaw. Once a 
    XMMEPIC instance is initialized for a given detector, filter, operation mode and date, 
    ECFs values can be obtained for different values of NH and photon index. Returned ECFs are 
    in units of counts × cm² / erg and it is possible to include the correction due
    to absorption.

    ECFs for the PN detector are very stable across epochs and operation modes, so it is 
    safe to use the default values. MOS detectors show higher variation across different 
    epochs, but still within ~2-3 per cent in the most extreme cases.

    The default values we use correspond to those assumed by [Rosen+2016](https://doi.org/10.1051/0004-6361/201526416) 
    for estimating fluxes for the XMM-Newton serendipitous catalogues.
    
    Examples
    --------
    - XMM ECF for PN in FullFrame with a Medium filter at energy band 3 (1-2 keV), assuming an
    absorbed powerlaw with a Hydrogen column density of 5×10²¹ cm-2 and photon index 1.9:
    
    >>> xmmpn_ecf = XMMEPIC("EPN", "Medium", eband="3")

    ECF with no absorption correction:

    >>> xmmpn_ecf(nh=5e21, gamma=1.9)
    <Quantity 5.61893163e+11 cm2 / erg>

    ECF including absorption correction:

    >>> xmmpn_ecf(5e21, 1.9, abscorr=True)
    <Quantity 3.27472149e+11 cm2 / erg>

    - ECF for the MOS2 camera with Thin filter at 2-10 keV for a recent observation, 
    using the default spectral parameters (NH = 3×10²⁰, Γ = 1.7):

    >>> xmmm2_ecf = XMMEPIC("EMOS2", "Thin", eband="HARD", date="2024-6-01")
    >>> xmmm2_ecf()
    <Quantity 4.41434307e+10 cm2 / erg>

    - Show available filters:
    
    >>> XMMEPIC.filters

    - Show available energy bands:
    
    >>> XMMEPIC.ebands

    - Show calibration epochs:
    
    >>> XMMEPIC.epochs
    """

    modes = {
        "epn": ("ff", "ef", "sw", "lw"),
        "emos": ("im",),
    }

    filters = ("Thin", "Medium", "Thick")

    ebands = {
        "1": (0.2, 0.5),
        "2": (0.5, 1.0),
        "3": (1.0, 2.0),
        "4": (2.0, 4.5),
        "5": (4.5, 12.0),
        "6": (0.2, 2.0),
        "7": (2.0, 12.0),
        "8": (0.2, 12.0),
        "9": (0.5, 4.5),
        "SOFT": (0.5, 2.0),
        "HARD": (2.0, 10.0),
    }

    # From https://www.cosmos.esa.int/web/xmm-newton/epic-response-files
    epochs = {
        "epn": {
            "e1": ("1999-12-10", "2007-01-01"),
            "e2": ("2007-01-01", "2014-01-01"),
            "e3": ("2014-01-01", "2021-01-01"),
            "e4": ("2021-01-01", Time(datetime.now())),
        },
        "emos": {
            "e1": ("1999-12-10", "2000-10-03"),
            "e2": ("2000-10-03", "2001-04-22"),
            "e3": ("2001-04-22", "2001-11-07"),
            "e4": ("2001-11-07", "2002-05-26"),
            "e5": ("2002-05-26", "2002-11-05"),
            "e6": ("2002-11-05", "2004-01-14"),
            "e7": ("2004-01-14", "2005-02-14"),
            "e8": ("2005-02-14", "2006-03-22"),
            "e9": ("2006-03-22", "2007-04-24"),
            "e10": ("2007-04-24", "2008-05-28"),
            "e11": ("2008-05-28", "2009-07-01"),
            "e12": ("2009-07-01", "2010-08-03"),
            "e13": ("2010-08-03", "2011-09-07"),
            "e14": ("2011-09-07", "2013-04-27"),
            "e15": ("2013-04-27", "2014-12-16"),
            "e16": ("2014-12-16", "2016-08-05"),
            "e17": ("2016-08-05", "2018-03-26"),
            "e18": ("2018-03-26", "2019-11-14"),
            "e19": ("2019-11-14", Time(datetime.now())),
        },
    }

    def __init__(self, detector, filter, eband="SOFT", mode=None, date=None):
        self._ecf = XMMECFValues()

        # Attributes set in the parse method
        self._parse_args(detector, filter, eband, mode, date)
        self._interpolators = self._set_interpolators()

    def _parse_args(self, detector, filter, eband, mode, date):
        self.detector = self._parse_detector(detector)
        self.filter = self._parse_filter(filter)
        self.mode = self._parse_mode(mode)
        self.eband = self._parse_eband(eband)
        self.epoch = self._parse_date(date)

    def _parse_detector(self, detector: str):
        try:
            return XMMDetector[detector]
        
        except KeyError:
            raise ValueError(f"Unknown detector: {detector}")

    def _parse_filter(self, filter: str) -> str:
        if filter == "Thin1" or filter == "Thin2" or filter == "Thin":
            filter = "Thin"

        if filter not in self.filters:
            raise ValueError(f"Unknown filter: {filter}")
        
        return filter
    
    def _parse_mode(self, mode: str) -> str:
        if mode is None:
            if self.detector.type == "epn":
                mode = "ff"
            else:
                mode = "im"

        if mode not in self.modes[self.detector.type]:
            raise ValueError(f"Unknown mode: {mode}")
        
        return mode
        
    def _parse_eband(self, eband: str) -> str:
        if eband not in self.ebands:
            raise ValueError(f"Unknown energy band: {eband}")
        
        return eband
        
    def _parse_date(self, date) -> str:
        if date is None:
            if self.detector.type == "epn":
                epoch = "e2"
            else:
                epoch = "e13"

        else:
            date = Time(date)
            epoch = None

            for key_epoch, values in self.epochs[self.detector.type].items():
                date_min = Time(values[0])
                date_max = Time(values[1])

                if date_min <= date <= date_max:
                    epoch = key_epoch
                    break

        if epoch is None:
            raise ValueError("Date is not compatible with the XMM mission.")
        
        return epoch

    def _set_interpolators(self):
        ecf_values_nocorr = np.array(
            self._ecf.nocorr[self.detector.tag][self.epoch][self.mode][self.eband][self.filter]
        )
        ecf_values_abscorr = np.array(
            self._ecf.abscorr[self.detector.tag][self.epoch][self.mode][self.eband][self.filter]
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


class XMMECFValues(metaclass=SingletonMeta):
    data_path = resources.files("ecfxa.data")

    def __init__(self) -> None:
        # Values with no Galactic absorption correction
        with resources.as_file(self.data_path / "xmm_ecfs.json.gz") as ecf_file:
            with gzip.open(ecf_file) as fp:
                self.nocorr = json.load(fp)

        # Values taking into account Galactic absorption correction
        with resources.as_file(self.data_path / "xmm_abscorr_ecfs.json.gz") as ecf_file:
            with gzip.open(ecf_file) as fp:
                self.abscorr = json.load(fp)


class XMMDetector(Enum):
    EMOS1 = ("mos1", "M1", "emos")
    EMOS2 = ("mos2", "M2", "emos")
    EPN = ("pn", "PN", "epn")

    def __init__(self, short, tag, type):
        self.short = short
        self.tag = tag
        self.type = type
