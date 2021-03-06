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

Features
--------

Here we gather usefull code that we daily use in the lab. 
TODO : complete

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
