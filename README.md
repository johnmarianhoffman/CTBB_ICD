# CTBB_ICD

## Are you sure you didn't mean FreeCT_ICD?

The "public" version of the software in this repository can be found at [FreeCT\_ICD](https://github.com/FreeCT/FreeCT_ICD).  This repository is for our "in-house" version and targets a reader library that is not publicly available.  I'm not saying you *can't* work from this project, but you may want the other one...

If you do want to work from here, you'll need to reimplement the reader functions and corresponding library.  There are several ways to do this, but the most portable is probably the following:

1. Change or reimplement: the "ctbb\_read.h" header file (included in setup.cpp)
2. Write: Your own ctbb\_read.cpp that reimplements the "void Read*Frame..." functions used to load raw data (these functions can be found in setup.cpp)
3. Compile: Compile your ctbb\_read.cpp into a ctbb\_read.a and put on the library load path
4. Modify the makefile: Link with -lctbb\_read (or whatever you change the name to)

## What is CTBB_ICD?

CTBB\_ICD is the model-based iterative reconstruction companion software to [CTBangBang](https://github.com/captnjohnny1618/CTBangBang).  It is free software (released under the GNU GPLv2) intended for use in research and education.  Modification and contributions are encouraged.

### Quick (Science) Facts:

CTBB_ICD...

* Reconstructs third-generation, helical, diagnostic, multidector CT data (without gantry tilt)
* Uses **iterative coordinate descent** as an optimizer
* Employs a **stored system matrix (SSM)** 
* Employs a **rotating slice** geometry to keep the stored system matrix sizes small (J. Xu and B. M. W. Tsui, “Iterative image reconstruction in helical cone-beam x-ray CT using a stored system matrix approach.,” Phys. Med. Biol., vol. 57, no. 11, pp. 3477–97, 2012.)
* Uses a **modified Joseph's method** approach for the projector
* Supports **flying focal spots** (*Note: at present, only an in-plane flying focal spot is supported, however support for Z and Z+in-plane are under active development*)

### Quick (Computing) Facts:

* Is **low-dependency**.  You only need to have the BOOST (https://www.boost.org/), yaml-cpp (https://github.com/jbeder/yaml-cpp), and [FreeCT_Reader](https://github.com/FreeCT/FreeCT_Reader) libraries installed on your system.
* Uses OpenMP for acceleration (customizable number of cores)
* Reconstructs a clinical CT dataset in about 12 hours

## What it is not

### CTBB_ICD is not what the manufacturers use

While our algorithms are based off of publications that are perhaps relevant to some of the current algorithms used in industry, they are not the algorithms used on clinical CT scanners and we make no claims to the similarity between our reconstructed images and what is arrived at clinically. In fact, moreso than FreeCT\_wFBP, FreeCT\_ICD employs a unique forward projector that has not been previously published.

Work has been done to objectively evaluate the quality of our reconstructions.  This can be found here:

*(insert link to technical note and/or whitepaper)* (**Coming soon**)

### CTBB_ICD is not a library

There are many great reconstruction libraries out there (http://conrad.stanford.edu/, http://www.openrtk.org/ to name two), and perhaps one day CTBB_ICD will be recast as a library.  Currently however, it is not a library, it a program.

It is structured modularly, so that major subsections of the reconstruction process are easy to identify/customize/edit/etc. so there are library-like qualities to the project to make it easy to use.

CTBB_ICD is designed to be compiled and run to reconstruct projection data from start to finish.

## Versions

The latest working version can be found on GitHub at https://github.com/captnjohnny1618/CTBB_ICD

Bleeding edge updates can be found at https://github.com/captnjohnny1618/CTBB_ICD/tree/develop

## Installation

*(to be added later)*

## Use

*(to be added later)*

## License

GNU GPLv2

*(more info to be added later)*

Copyright 2018 John Hoffman, Stefano Young, Frédéric Noo
