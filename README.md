# Heliumtools 
![Python3.9](https://img.shields.io/badge/python-3.9-blue)
![Python3.10](https://img.shields.io/badge/python-3.10-red)


This package is developped by the [Quantum Atom Optics]((https://www.lcf.institutoptique.fr/groupes-de-recherche/gaz-quantiques/experiences/quantum-atom-optics)) group at the  Laboratoire Charles Fabry, [Institut d'Optique Graduate School](https://www.institutoptique.fr/) in Palaiseau (France). 

We make use of Bose-Einstein condensates (BEC) of metastable helium atom (He*) to perform interferometry experiments inspired by quantum optics. A remarkable feature of the metastable state of helium is its very large internal energy (20eV), allowing the use of a single-atom resolved detector, based on the electron amplification. This detector is called a microchannel plate (MCP) and provides three-dimensional data in the momentum space of the atoms. This python package gathers code we daily use to analyze our datas. 



- Free software: MIT license

<a href = "https://github.com/Tanu-N-Prabhu/Python/graphs/contributors">
  <img src = "https://contrib.rocks/image?repo = GitHub_username/repository_name"/>
</a>

Made with [contributors-img](https://contrib.rocks).

## Installation

This package was written in Python 3.9, so we recommend using it with Python 3.9 or newer versions.

- Create a new environment (heliumenv here for example). We use mkvirtualenv.

    - `mkvirtualenv -p /usr/bin/python3.9 heliumenv`
    - If it is not done, add this shortcut to your path bashrc: `gedit ~/.bashrc` and paste:

        ```bash
        export WORKON_HOME=$HOME/.virtualenvs
        export PROJECT_HOME=$HOME/Devel
        source /usr/share/virtualenvwrapper/virtualenvwrapper.sh
        ```

- Clone the heliumtools repository and pip install it `cd /your/path/heliumtools && pip install -e .`
- Since March 2023, we use PyTorch to compute correlations. We add trouble to install it through pip but one can install it using `pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu`

## Modules
### Momentum resolved correlations : the correlations module


### Momentum integrated correlations : the correlations2 module

This module ables one to probe integrated correlations. See [the module documentation](manual_doc/correlations2.md) for more informations.
<div style="text-align: center;">
    <img src="manual_doc/img/correlations2_local.png" alt="Local correlations." width="600">
</div>


### ODT-calc : trap properties calculations

This module analyze cold-atom traps combining (far off-resonant) laser beams and (circular) magnetic coils. It was first developped by [Alexandre Dareau](https://github.com/adareau) in the [odt-calc package](https://github.com/adareau/odt-calc) but was forked and implemented into the heliumtools package in 2022. 


<div style="text-align: center;">
    <img src="manual_doc/img/odt_calc0.png" alt="ODT-calc image should appear here." width="600">
</div>

*This image shows the analysis that can be carried by the odt-calc module. It is possible to model laser beams and magnetic coils.*


<!-- Here we gather useful code that we use daily in the lab. 
One will recover:
- Correlations: is the code we use to check whether or not atomic pairs are correlated,
- Odt-calc: is the code we use to check properties of our optical dipole trap or magnetic trap,
- Bragg: this class enables the user to compute reflectivity profiles to set up our Bragg pulses. -->

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [`audreyr/cookiecutter-pypackage`](https://github.com/audreyr/cookiecutter-pypackage) project template.
