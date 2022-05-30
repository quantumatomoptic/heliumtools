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
To use all fonctioncalities of this package, you will need to install HAL from https://github.com/adareau/HAL. To do so:

- create a new environnement (heliumenv here for exemple). We use mkvirtualenv.

    - ``mkvirtualenv -p /usr/bin/python3.9 heliumenv``
    - if it is not done, add this shortcut to your path bashrc  : ``gedit ~/.bashrc``   and paste::
    
        export WORKON_HOME=$HOME/.virtualenvs
        export PROJECT_HOME=$HOME/Devel
        source /usr/share/virtualenvwrapper/virtualenvwrapper.sh
- install PyQt5 and HAL in the environnent. 

    - ``pip install PyQt5`` : we must install PyQt5 before HAL to prevent issues.
    - ``cd ~/HAL && pip install -e .`` sot that you will install HAL as a developper. 
    
- clone the heliumtools repository and pip install it ``cd /your/path/heliumtools && pip install -e .``

Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
