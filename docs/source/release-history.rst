===============
Release History
===============


Release V0.3.0  (2022-02-23)
----------------------------

    Applications:
        + Cathode Viewer
        + FAR Viewer
        + Image Viewer

    cathode module:
        + Methods to laod Cathode Camera images and perform aspect corrections.
        + Methods to calculate temperature based on emissivity and camera QE.
        + Method to calculate cathode temperature based on power model of current and voltage supplied.

    far_pyro module:
        + Methods to load FAR pyrometer data
        + Method to calculate the approximate temperature from a two-band parameterized technique
        + Planck temperature fit to pyrometer data with calibration and emissivity functions

    planck module:
        + Planck distribution
        + Planck spectral radiance pdf distribution and other Planck fit methods
        + Powerlaw distribution and fit methods
        + Stefan Boltzman equation

    Other modules:
        + Image loading
        + Controls information methods
        + Misc Utilities
