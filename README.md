TurbospectrumPy
===

`This repository is maintained by` **J. G. Fernandez-Trincado**. `Feel free to check it out, make comments, and provide me with some feedback. I will update the design of it eventually, since at the moment, I have just literally slapped a beta version onto the "default" design of the repository.`

If you have any comments, suggestions, and/or complaints, please contact me: jfernandez@obs-besancon.fr and/or jfernandezt87@gmail.com

This repository is a simple Python Script to run TurboSpectrum
--

  * Turbospectrum code for spectral synthesis v15.1 can be found at [Turbospectrum code for spectral synthesis v15.1](http://www.pages-perso-bertrand-plez.univ-montp2.fr)
  * I also recommend to use a more complete version for spectral synthesis for [Apogee](https://github.com/Fernandez-Trincado/apogee) spectra, developed in Python by Jo Bovy.

Manual analysis 
--

[Manualsynthesis.py](https://github.com/Fernandez-Trincado/TurbospectrumPy/blob/master/Manualsynthesis.py) is a Python script for the treatment and analysis of stellar spectra. Some of the main functionalities of this script are the following:

![Figure 1](https://github.com/Fernandez-Trincado/TurbospectrumPy/blob/master/Abundances2.png)


Instructions to make spectrum synthesis star-to-star and line-to-line. Below are instructions for this:
--

 0. NOTE: This program only works for APOGEE spectra, but can be extended and other cases.

 0. Download the main Python program at: [Manualsynthesis.py](https://github.com/Fernandez-Trincado/TurbospectrumPy/blob/master/Manualsynthesis.py)

 1. Download the atmosphere model at:
 
   1.1 [Kurucz model](https://data.sdss.org/sas/apogeework/apogee/spectro/redux/speclib/kurucz_filled/mm08cp00op00/) (this requires your trac account) 
   
   1.2 [MARCS model](https://data.sdss.org/sas/apogeework/apogee/spectro/redux/speclib/marcs/edvarsson/) (this requires your trac account)
 
 2. Download your spectrum at:
  
   2.1 [SDSS-IV Spectra](https://data.sdss.org/sas/apogeework/apogee/spectro/redux/r6/stars/l30e/l30e.2/) (this requires your trac account)
 
 3. Run the Python script:
 
  3.1 Option 1: > python Manualsynthesis.py
  
  3.2 Option 2: > ./Manualsynthesis.py
  
 4. Input and Output (example). This program find the best-fit between the synthetic spectrum and the observed spectrum using the Chi^2 technique.
 
    ![Figure 2](https://github.com/Fernandez-Trincado/TurbospectrumPy/blob/master/run.png)
 
 5. `Disclaimer:` The manual analysis is made at your own risks. The author is not responsible for wrong applications of their program. In case you would like complementary informations, you may contact me directly (jfernandez[at]obs-besancon.fr and/or jfernandezt87[at]gmail.com). 

Instructions to make spectrum synthesis for many APOGEE stars line-to-line. Below are instructions for this:
--

In construction ...


References
--

 1. Conversion from vacuum to standard air wavelengths can be found at [VAC2AIR/AIR2VAC](http://hebe.as.utexas.edu/apogee/docs/air_vacuum.pdf)
 2. Radiative Transfer in Stellar Atmospheres can be found at [Review](http://www.staff.science.uu.nl/~rutte101/rrweb/rjr-edu/coursenotes/rutten_rtsa_notes_2003.pdf) 
  
Support TurbospectrumPy
--

If you use TurbospectrumPy in your research, we would be grateful if you could include an acknowledgment in papers and/or presentations:

    "This research made use of TurbospectrumPy, a community-developed core Python package for Astronomy."

Thanks.

  
  
  
  
