nir-spectrum-analysis
=====================

Welcome to the NIR spectrum analysis software package, written for MATLAB.

CREDIT: The intitial code base was created by Qiong Zhang as part of his M.Sc. in Dr. Jeff Dunn's lab. Further modifications were done by Thomas Johnson (main routines) and Daniel Korchinski (cytochrome oxidase spectrum).

Instructions:
1. Turn on the andor spectrometer and ensure it is functioning b using the andor software package.
2. Add this folder and its subfolders to your MATLAB path.
3. Run 'main.m'. This is the gui initialization routine that will then execute 'part.m'.
4. Initialize the system ('Initialize' button), setting the target temperature to -40 degrees Celcius. This ensures that the CCD chip is adequately cooled to reduce dark noise and calm down any bad pixels.
5. Spectral calibration: this is done using the known spectrum of a neon lamp. One is located on the Oriel light source (the little orange light). After measuring the spectrum, you must selectt he peaks of the measured spectrum (red) that correspond to the theoretical spectrum (blue). You click on the red peak, then its corresponding blue peak. After all the peaks have been done, press enter to finish.
6. Intensity calibration: in order to determine the attenuation that occurs at each wavelength, we must first learn what the original light intensity going into the tissue of interest is at each wavelength. This is done by taking a measurement of the light coming out of the transmitting fibre with the receiving fibre. Due to the strength of this signal, the optical density filters must be used to attenuate the signal enough to avoid the CCD chip being saturated. A setting of 5 (meaning a 10^5 reduction in light intensity by the filter) is usually sufficient.
7. Measurement! Follow the on-screen instructions and use the following default settings for good results:
- number of averages per measurement: 10
- length of time per measurement: 0.5 seconds



