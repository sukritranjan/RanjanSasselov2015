This folder contains a Python implementation of the formalism described in Section 3.4 and Appendix C of Ranjan & Sasselov (2015, Astrobiology, accepted). If you make use of this code, please cite this work. 

Comments, criticisms and suggestions welcomed. Please email feedback to: sranjan@cfa.harvard.edu.

PURPOSE
The formalism implemented here estimates the column density of CO2 required to quench photolysis to the degree required to permit HCN and CH4 levels to build up to a given level. 

These gases are speculated to have been important feedstock gases for prebiotic synthesis on the early Earth. However, few empirical limits are available on their abundance. Both gases are strong absorbers in the UV, and photolysis is expected to be a key constraint on their abundances. Here, we implement a simple model to estimate the surface partial pressures of these gases assuming a single source (geochemical flux) and a single sink (photolysis).

For purposes of this code, we assume a 1-bar, 290 K isothermal atmosphere in hydrostatic equilibrium. We assume the atmosphere to be CO2/N2 dominated, and for all gases to be well-mixed until a height z0 at which the feedstock gas mixing ratio goes to 0. We ignore self-shielding; hence, CO2 column densities estimated here are upper bounds. 

We emphasize that these calculations are *not* intended as realistic models of the early Earth's atmosphere. Rather, our intention is to enable those interested in prebiotic chemistry to constrain to first order the levels of CH4 and HCN available for synthesis reactions as a function of assumed CO2 partial pressure. 



FILE DESCRIPTIONS
---PhotolysisCalculation.py: Core implementation of Ranjan & Sasselov (2015) formalism used to estimate level of CO2 required to shield buildup of CH4 and HCN to a given surface pressure, in the form of functions that estimate the CH4/HCN partial pressure as a function of N_CO2 and z0.

---CO2_HuestisBerkowitz(2010)_300K_0.1254-201.6nm(evaluation).txt: Absorption cross-sections of CO2. Taken from http://satellite.mpic.de/spectral_atlas/cross_sections/Carbon-oxides/CO2_HuestisBerkowitz%282010%29_300K_0.1254-201.6nm%28evaluation%29.txt; based on a compilation originally reported in Huestis & Berkowitz (2010)

---CH4_Au(1993)_298K_5.6-165nm(e,e).txt: Absorption cross-sections of CH4. Taken from http://satellite.mpic.de/spectral_atlas/cross_sections/Alkanes+alkyl%20radicals/Alkanes/CH4_Au%281993%29_298K_5.6-165nm%28e,e%29.txt; based on measurements originally reported in Au et al (1993)

---HCN_cxs.dat: Absorption cross-sections of HCN. Originally reported by Nuth et al (1982). This electronic version was kindly compiled and made available to us by V. Vuitton (1/12/2015)

---claire_3d9Ga.dat: Incident flux from the 3.9 Ga sun at 1 AU. Generated from the code of Claire et al (2012), which is available at: http://depts.washington.edu/naivpl/content/models/solarflux

---cookbook.py: Compendium of python functions useful for manipulations involving binned data. Contains both functions from scipy cookbooks as well as homebrew code.