===========
heliumtools
===========


.. image:: https://img.shields.io/pypi/v/heliumtools.svg
        :target: https://pypi.python.org/pypi/heliumtools

.. image:: https://img.shields.io/travis/quantumatomoptic/heliumtools.svg
        :target: https://travis-ci.com/quantumatomoptic/heliumtools

.. image:: https://readthedocs.org/projects/heliumtools/badge/?version=latest
        :target: https://heliumtools.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Some tools for the Helium1 team @LCF. We are a group of physicist working on ultracold Helium at Laboratoire Charles Fabry in Palaiseau (France). Check out our webpage: https://www.lcf.institutoptique.fr/groupes-de-recherche/gaz-quantiques/experiences/quantum-atom-optics


* Free software: MIT license
* Documentation: https://heliumtools.readthedocs.io.


Installation
--------
This package was written in Python3.9 so we recommand to use it with python3.9. On Linux, one can 

- create a new environnement (heliumenv here for exemple). We use mkvirtualenv.

    - ``mkvirtualenv -p /usr/bin/python3.9 heliumenv``
    - if it is not done, add this shortcut to your path bashrc  : ``gedit ~/.bashrc``   and paste::
    
        export WORKON_HOME=$HOME/.virtualenvs
        export PROJECT_HOME=$HOME/Devel
        source /usr/share/virtualenvwrapper/virtualenvwrapper.sh
- clone the heliumtools repository and pip install it ``cd /your/path/heliumtools && pip install -e .``
- since March 2023, we use PyTorch to compute correlations. We add trouble to install it through pip but one can install it using ``pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu``

Features
--------

Here we gather usefull code that we daily use in the lab. 
One will recover : 
* correlations : is the code we use to check wether or not atomic pairs are correlated,
* odt-calc : is the code we use to check properties of our optical dipole trap or magnetic trap,
* Bragg : this class ables the user to compute reflexivity profiles to setup our Bragg pulses.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
