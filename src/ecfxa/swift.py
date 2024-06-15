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


class SwiftXRT:
    """
    Energy Conversion Factors (ECFs) for the XRT instrument on-board the Swift telescope.

    We include ECF estimates for the two different operation modes of XRT (wt and pc),
    their corresponding grades (0 and 0-2 for wt, 0, 0-4 and 0-12 for pc), and the different
    epochs defined by the Swift calibration team. ECFs were estimated assuming an absorbed 
    powerlaw. Once a SwiftXRT instance is initialized for a given mode, grade, energy band 
    and date, ECFs values can be obtained for different values of NH and photon index. 
    Returned ECFs are in units of counts × cm² / erg and it is possible to include the 
    correction due to absorption.
    
    Examples
    --------
    - XRT ECF for pc mode with event grade 0 to 4 at energy band 2 (1-2 keV), assuming an
    absorbed powerlaw with a Hydrogen column density of 5×10²¹ cm-2 and photon index 1.9:

    >>> xrtpc_ecf = SwiftXRT("pc", grade="04", eband="2")

    ECF with no absorption correction:

    >>> xrtpc_ecf(nh=5e21, gamma=1.9)
    <Quantity 5.20816877e+10 cm2 / erg>
    
    ECF including absorption correction:

    >>> xrtpc_ecf(5e21, 1.9, abscorr=True)
    <Quantity 3.03442833e+10 cm2 / erg>
    
    - Show grades for the different modes:
    
    >>> SwiftXRT.grades
    
    - Show available energy bands:
    
    >>> SwiftXRT.ebands

    - Show calibration epochs:

    >>> SwiftXRT.epochs
    """

    grades = {
        "wt": ("0", "02"),
        "pc": ("0", "04", "012"),
    }

    ebands = {
        "0": (0.3, 10.0),
        "1": (0.3, 1.0),
        "2": (1.0, 2.0),
        "3": (2.0, 10.0),
        "SOFT": (0.5, 2.0),
        "HARD": (2.0, 10.0),
    }

    # From https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/swift/docs/xrt/SWIFT-XRT-CALDB-09_v22.pdf
    epochs = {
        "e1": ("2004-12-01", "2007-01-01"),
        "e2": ("2007-01-01", "2007-08-31"),
        "e3": ("2007-08-31", "2009-01-01"),
        "e4": ("2009-01-01", "2011-01-01"),
        "e5": ("2011-01-01", "2013-01-01"),
        "e6": ("2013-01-01", "2013-12-12"),
        "e7": ("2013-12-12", "2021-01-01"),
        "e8": ("2021-01-01", Time(datetime.now())),
    }

    def __init__(self, mode, grade="0", eband="SOFT", date=None):
        self._ecf = SWXRTECFValues()

        # Attributes set in the parse method
        self._parse_args(mode, grade, eband, date)
        self._interpolators = self._set_interpolators()

    def _parse_args(self, mode, grade, eband, date):
        self.mode = self._parse_mode(mode)
        self.grade = self._parse_grade(grade)
        self.eband = self._parse_eband(eband)
        self.epoch = self._parse_date(date)


    def _parse_mode(self, mode: str):
        if mode not in self.grades:
            raise ValueError(f"Unknown mode: {mode}")

        return mode            
    
    def _parse_grade(self, grade: str) -> str:
        if grade not in self.grades[self.mode]:
            raise ValueError(f"Unknown grade: {grade}")
        
        return grade
        
    def _parse_eband(self, eband: str) -> str:
        if eband not in self.ebands:
            raise ValueError(f"Unknown energy band: {eband}")
        
        return eband
        
    def _parse_date(self, date) -> str:
        if date is None:
            epoch = "e6"

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
            raise ValueError("Date is not compatible with the Swift mission.")
        
        return epoch

    def _set_interpolators(self):
        ecf_values_nocorr = np.array(
            self._ecf.nocorr[self.mode][self.epoch][self.grade][self.eband]
        )
        ecf_values_abscorr = np.array(
            self._ecf.abscorr[self.mode][self.epoch][self.grade][self.eband]
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


class SWXRTECFValues(metaclass=SingletonMeta):
    data_path = resources.files("ecfxa.data")

    def __init__(self) -> None:
        # Values with no Galactic absorption correction
        with resources.as_file(self.data_path / "swift_ecfs.json.gz") as ecf_file:
            with gzip.open(ecf_file) as fp:
                self.nocorr = json.load(fp)

        # Values taking into account Galactic absorption correction
        with resources.as_file(self.data_path / "swift_abscorr_ecfs.json.gz") as ecf_file:
            with gzip.open(ecf_file) as fp:
                self.abscorr = json.load(fp)
