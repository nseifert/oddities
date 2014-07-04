Splatalogue Auto-Assigner 
==========================

This is a simple script written in Python 2.7x that, given an input Fourier transformed spectrum 
[with columnal data (freq, inten)], will peakpick the spectrum within 
a given intensity threshold and return the N closest Splatalogue results
for each peakpicked frequency. 

There is additional information under __main__ on usage instructions.

Library prerequisties:
* Numpy 
* Requests
* BeautifulSoup
