# OpenBiology_2016

This repository, ideally, contains the necessary code and data to re-create the results of Hockenberry et al. 2017 ( see: http://rsob.royalsocietypublishing.org/content/7/1/160239.long ). Some pieces of data, however, may be excessively large to include here under which circumstances they may be ommitted and the users should reach out via the contact information listed below. 

The code is primarily in the form of Jupyter Notebooks to:

a) create translation efficiency measurements from ribosome profiling and RNA-seq .wiggle data (provided in Data/) (Please note that we did not perform any sequence mapping on this data and instead relied on the original authors mappings and .wig files as our intention was to compare as closely as possible to their results. Ideally, this mapping should probably be re-run especially in the event of course that the user is hoping to test their own data and given newer mapping / analysis protocols that are better able to handle the specifics of ribosomal profiling data). Users should first run through `calculate_trans_eff.ipynb' if they want to see how we processed and worked with ribosomal profiling and mRNA seq data. The major intention of this notebook however is to create .json files of translation efficiency measurements and save them to ../Data/. We currently provide these files to that the user does not have to re-run this notebook at all should they choose not to.   

b) to use that translation efficiency data in order to investigate the effect of anti-Shine-Dalgarno sequence binding. A variety of files are included in Data/ that quantify aspects of cis-mRNA folding and trans-anti-Shine-Dalgarno sequence binding strengths.

The code should run with a basic Anaconda running Python 3.xx or similar install that includes basic scientific libraries. Users can check the top of .py libraries or .ipynb jupyter notebooks for individual requirements but basic scipy,numpy,pandas,etc. shoud do the trick. The exception being that Biopython should be installed which I leave as an exercise to the user.

Specific questions/comments/concerns should be addressed to Adam Hockenberry via adam [dot] hockenberry [at] northwestern [dot] edu. I'll try and get back to you in a timely manner.
