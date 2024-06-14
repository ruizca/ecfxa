import os
from datetime import datetime
from pathlib import Path

import pxsas


def make_ccf(ccf_path, date=None):
    if date is None:
        date = datetime.now()
        date = date.isoformat()
    
    pxsas.run(
        "cifbuild",
        calindexset=ccf_path,
        withobservationdate=True,
        observationdate=date,
        # analysisdate=date,
    )


data_path = Path("data", "xmm", "tests")

kwargs_arfgen_onaxis = {
    "withrmfset": False,
    "withsourcepos": True,  # on-axis position (no vignetting correction)
    "sourcecoords": "tel",
    "sourcex": 0,
    "sourcey": 0,
    "filterdss": False,
    "withdetbounds": True,
    "withbadpixcorr": False,  # no correction of bad pixels
    "modelee": False,    # no correction of PSF enclosed energy fraction
    "detmaptype": "flat", 
    "detxbins": 1,
    "detybins": 1,
    # "applyxcaladjustment": True,
    # "applyabsfluxcorr": True,
}


# Generate ARFs for pn, mos1 and mos2 for two different spectra, assuming that the observation was done today
ccf_path = data_path / "nowccf.cif"
# make_ccf(ccf_path)
os.environ["SAS_CCF"] = str(ccf_path.resolve())

# Following http://xmm-tools.cosmos.esa.int/external/sas/current/doc/arfgen/node11.html
pxsas.run(
    "arfgen",
    arfset=data_path / "now_PNS13.arf",
    spectrumset=data_path / "PNS13SRSPEChz403409.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_PNS16.arf",
    spectrumset=data_path / "PNS16SRSPEChz421525.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_M1S13.arf",
    spectrumset=data_path / "M1S13SRSPEChz403409.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_M1S16.arf",
    spectrumset=data_path / "M1S16SRSPEChz421525.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_M2S13.arf",
    spectrumset=data_path / "M2S13SRSPEChz403409.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_M2S16.arf",
    spectrumset=data_path / "M2S16SRSPEChz421525.FTZ",
    **kwargs_arfgen_onaxis,
)


pxsas.run(
    "arfgen",
    arfset=data_path / "now_PNS13_medium.arf",
    spectrumset=data_path / "PNS13SRSPEChz403409_medium.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_PNS13_thick.arf",
    spectrumset=data_path / "PNS13SRSPEChz403409_thick.FTZ",
    **kwargs_arfgen_onaxis,
)

pxsas.run(
    "arfgen",
    arfset=data_path / "now_M1S13_medium.arf",
    spectrumset=data_path / "M1S13SRSPEChz403409_medium.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_M1S13_thick.arf",
    spectrumset=data_path / "M1S13SRSPEChz403409_thick.FTZ",
    **kwargs_arfgen_onaxis,
)

pxsas.run(
    "arfgen",
    arfset=data_path / "now_M2S13_medium.arf",
    spectrumset=data_path / "M2S13SRSPEChz403409_medium.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_M2S13_thick.arf",
    spectrumset=data_path / "M2S13SRSPEChz403409_thick.FTZ",
    **kwargs_arfgen_onaxis,
)


pxsas.run(
    "arfgen",
    arfset=data_path / "now_PNS13_singles.arf",
    spectrumset=data_path / "PNS13SRSPEChz403409_singles.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_PNS13_medium_singles.arf",
    spectrumset=data_path / "PNS13SRSPEChz403409_medium_singles.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "now_PNS13_thick_singles.arf",
    spectrumset=data_path / "PNS13SRSPEChz403409_thick_singles.FTZ",
    **kwargs_arfgen_onaxis,
)


# Generate ARFs for pn, mos1 and mos2 for two different spectra, assuming that the observation was done in 2001
ccf_path = data_path / "2001_ccf.cif"
# make_ccf(ccf_path, date="2001-01-01")
os.environ["SAS_CCF"] = str(ccf_path.resolve())

# Following http://xmm-tools.cosmos.esa.int/external/sas/current/doc/arfgen/node11.html
pxsas.run(
    "arfgen",
    arfset=data_path / "2001_PNS13.arf",
    spectrumset=data_path / "PNS13SRSPEChz403409.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "2001_PNS16.arf",
    spectrumset=data_path / "PNS16SRSPEChz421525.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "2001_M1S13.arf",
    spectrumset=data_path / "M1S13SRSPEChz403409.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "2001_M1S16.arf",
    spectrumset=data_path / "M1S16SRSPEChz421525.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "2001_M2S13.arf",
    spectrumset=data_path / "M2S13SRSPEChz403409.FTZ",
    **kwargs_arfgen_onaxis,
)
pxsas.run(
    "arfgen",
    arfset=data_path / "2001_M2S16.arf",
    spectrumset=data_path / "M2S16SRSPEChz421525.FTZ",
    **kwargs_arfgen_onaxis,
)
