# DDrates

[![arXiv](https://img.shields.io/badge/arXiv-2006.02453-B31B1B.svg)](https://arxiv.org/abs/2006.02453)
[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)


This is a python code for calculating the rates for nuclear recoil, Migdal effect and photon Bremsstrahlung for argon detectors using the Standard Halo Model. This code was used to compute the DM rates in [arXiv:2006.02453](https://arxiv.org/abs/2006.02453). The statistical analysis was performed with the code in https://github.com/piacent/LAr-MigdalLimits/blob/master/LArAnalysis.R.ipynb.



### Features

- Spin-independent dark matter - nucleus scattering for argon detectors, in particular DarkSide-50;
- Migdal effect and photon Bremsstrahlung for argon detectors;

### Contents

- `data/bremsstrahlung` - Contains the data for the X-ray form factors, from https://physics.nist.gov/PhysRefData/FFast/html/form.html;
- `data/DarkSide50` - Contains the data for the calibration curve to convert electron or nuclear recoil spectra to ionization spectra in LAr underground detectors, from 1802.06994 and 1802.06998;
- `data/Migdal` - Contains the differential transition probabilities for argon, from 1707.07258;
- `src` - Python files to compute the nuclear recoil, Migdal and Bremsstrahlung rates, plus an example jupyter notebook.

### Requirements

- numpy
- scipy
- numericalunits
- pandas

### Citing DMrates

If you use the code, please cite [arXiv:2006.02453](https://arxiv.org/abs/2006.02453):

>Grilli di Cortona G.,  Messina A., and Piacentini S. "Migdal effect and photon Bremsstrahlung: improving the sensitivity to light dark matter of liquid argon experiments", arXiv:2006.02453

as well as all the relevant articles for each module, listed in the files References.txt in the data/<name> folder.

### Authors

The author of this code is Giovanni Grilli di Cortona.

For questions, bug reports or suggestions please contact [ggrillidc@fuw.edu.pl](mailto:ggrillidc@fuw.edu.pl)
