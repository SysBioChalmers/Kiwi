
===========
Quick start
===========
For a proper installation of *Kiwi*, follow the instructions in the Package Documentation.

For a quick flavor of *Kiwi*, you can download the package from the Python Package Index (PyPI), open the command-prompt, navigate to the folder where the package was downloaded, and run:

	> tar -xvf KiwiDist-put.version.number.tar.gz

	> cd  KiwiDist-put.version.number

	> python -c "import kiwi; kiwi.plot(gsn='kiwi/demo/GSN.txt', gss='kiwi/demo/GSS_LUAD.txt', gls='kiwi/demo/GLS_LUAD.txt', gsc='kiwi/demo/GSC_LUAD.txt', nwf='exampleNetworkPlot.pdf', hmf='exampleHeatmap.pdf')"

This command fails if the dependencies (chiefly numpy, matplotlib, networkx) are not properly installed for the Python interpreter called by the command prompt. If you are using Canopy, follow the instructions above and precede the last command with: %system.

This run will produce two example figures. We recommend to follow proper installation procedures and use the Python interpreter to run *Kiwi*.