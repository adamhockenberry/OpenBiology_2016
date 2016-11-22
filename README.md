# OpenBiology_2016

This repository, ideally, contains the necessary code and data to re-create the results of Hockenberry et al. xxxx. pending acceptance of this manuscript.

The code is actively being compiled at the moment (as of Nov. 21 2016) so may not fully work for all users as some further data files remain to be added in the coming days/weeks.

The code is primarily in the form of Jupyter Notebooks to a) create translation efficiency measurements from ribosome profiling and RNA-seq .wiggle data (provided in Data/) and b) to use that translation efficiency data in order to investigate the effect of anti-Shine-Dalgarno sequence binding. A variety of files are included in Data/ that quantify aspects of cis-mRNA folding and trans-anti-Shine-Dalgarno sequence binding strengths.

The code should run with a basic Anaconda running Python 3.xx or similar install that includes basic scientific libraries. Users can check the top of .py libraries for individual requirements but basic scipy,numpy,pandas,etc. shoud do the trick. The exception being that Biopython should be installed which I leave as an exercise to the user.

More comments forthcoming as this code cleaning/compilation nears completion.
