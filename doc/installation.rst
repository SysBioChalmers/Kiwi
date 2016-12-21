Installation
================================

*Kiwi* can be found at the `Python Package Index (PyPI) <https://pypi.python.org/pypi/KiwiDist>`_ .

Requirements
------------

*Kiwi* was developed in Python 2.7 but should work also with Python 3. Other than Python, it depends on the following packages:

* matplotlib >= 1.3.1, <=1.4.3
* mygene >= 2.1.0
* networkx >= 1.8.1
* numpy >= 1.8.0
* pandas >= 0.13.1
* scipy >= 0.13.3, <=0.16.0
* six >=1.5

Package installation
--------------------
There a three ways to install packages from PyPI. More information can be found in `PyPI <https://wiki.python.org/moin/CheeseShopTutorial>`_.

Pip
^^^
This way is recommended to install dependencies automatically, if the Python path is set correctly in your system. Get to the command prompt and type: ::

	> pip install KiwiDist

Distutils
^^^^^^^^^
This package comes in bundle with Python. `Download <https://pypi.python.org/pypi/KiwiDist>`_ the package, extract it, get to the command prompt, navigate to the root directory, and type: ::

	> python setup.py install

The distribution will be installed into `site-packages` directory of the Python interpreter used to run the `setup.py install` command. No dependencies will be installed! Also, make sure to have admin privilege for installation, otherwise act as root user: ::

	> sudo python setup.py install

To check if the installation worked correctly, open the Python interpreter and type: ::

	> import kiwi

If no error is raised, then you are good to go!

To get started, open the Python interpreter and type: ::

        > import kiwi
        > ?kiwi.plot

Source
^^^^^^
For a quick flavor of *Kiwi*, you can download the package from the `Python Package Index (PyPI) <https://pypi.python.org/pypi/KiwiDist>`_,
open the command-prompt, navigate to the folder where the package was downloaded, and run: ::

	> tar -xvf KiwiDist-put.version.number.tar.gz
	> cd  KiwiDist-put.version.number
	> python -c "import kiwi; kiwi.plot(gsn='kiwi/demo/GSN.txt', gss='kiwi/demo/GSS_LUAD.txt', gls='kiwi/demo/GLS_LUAD.txt', gsc='kiwi/demo/GSC_LUAD.txt', nwf='exampleNetworkPlot.pdf', hmf='exampleHeatmap.pdf')"

This command fails if the dependencies (chiefly numpy, matplotlib, networkx) are not properly installed for the Python interpreter called by the command prompt. If you are using Canopy, follow the instructions above and precede the last command with: %system.

Package update
--------------
There are two recommended ways to check for updates, depending on which environment was used to install *Kiwi*.

Pip
^^^
Open the command prompt and simply issue: ::

	> pip install KiwiDist --upgrade

Canopy
^^^^^^
If Kiwi was first installed in `Canopy <https://www.enthought.com/products/canopy/>`_, to update the library in this environment, issue the following commands in the prompt window: ::

	%system pip install KiwiDist --upgrade (-U --no-deps, if no dependencies should be upgrade)

To check the current version: ::

	%system pip show KiwiDist

To check the current path: ::

	import kiwi; kiwi.__path__

.. image:: kiwi_logo.png
	:width: 200
	:align: right
