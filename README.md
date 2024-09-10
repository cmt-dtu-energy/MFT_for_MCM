# Mean field theory (MFT) model for magnetocalorics

This repository contains an implementation in Matlab of the mean field theory (MFT) model, also called the mean field model (MFM), for generating material data for magnetocaloric materials (MCMs).

The program takes as input the Curie temperature, the Landé factor, the total angular momentum, the Debye temperature, the number of spins pr unit mass, the molar mass and the Sommerfeld constant.
For a specified temperature range and magnetic field range, the program outputs the specific heat, the entropy, the magnetization and the adiabatic temperature change of the material.


## References and theory
The theory of the model can be found on p.6 (sec. 2.2.1) in the thesis found <A HREF="https://orbit.dtu.dk/files/101639258/Designing_a_magnet.pdf">here</A> or in [Smith, A., Bahl, C.R.H., Bjørk, R., Engelbrecht, K., Nielsen, K. K. and Pryds, N., "Materials challenges for high performance magnetocaloric refrigeration devices." Advanced Energy Materials 2, no. 11 (2012): 1288-1318, <A HREF="https://doi.org/10.1002/aenm.201200167">DOI: 10.1002/aenm.201200167</A>]. Please cite these works as well as this repository when utilizing the model to generate data.

## Acknowledgement
This work acknowledge funding from the EU-funded project HyLICAL, which received funding from the Europeans Union’s Horizon Europe research and innovation program under grant agreement No 101101461. The project is supported by the Clean Hydrogen Partnership and its members. We also acknowledge funding from the Norwegian Research Council project 336403 entitled "Hydrogen Liquefaction with Caloric materials" (LIQUID-H).
