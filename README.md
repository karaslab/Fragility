# Fragility
Fragility module for RAVE

This is an implementation of Li et al's [neural fragility](https://doi.org/10.1038/s41593-021-00901-w) biomarker of the epileptogenic zone within the RAVE environment. The eventual goal is to make computationally complex iEEG biomarkers of the epileptogenic zone such as neural fragility more accessible to epilepsy surgery programs and researchers.

RAVE is a powerful, free, open-source, NIH-funded research tool designed to analyze and visualize iEEG data that allows users to code their own modules for custom analyses. All information about RAVE can be found on [the RAVE website](https://rave.wiki).

## Installation Instructions

The Fragility module can be installed with the following command: `remotes::install_github('karaslab/fragility')`

After installing, it needs to be recognized by RAVE. Run `rave::rave_options()` and then click the "check for new modules button". The fragility module should now appear in the module list. Then, run `rave::arrange_modules` and restart R and RStudio. The next time you run `start_rave()`, the Fragility module will show up in the left-hand menu.
