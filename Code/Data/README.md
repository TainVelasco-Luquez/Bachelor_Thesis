Description of  the R objects inside each file linked to the HIassembly .Rmd, HIvalidation .Rmd and RBDmodule.Rmd, in order to make the rendering of such .Rmd faster:


0. RBDmodule.RData
Contains the most computationally expensive objects created in the file RBDmodule.RMD:
	0.0. igraph and network RBDmodule objects
	0.1. ERGM fit
	0.2. NEAT
	0.3. Transitivity and Averege Path Length (APL) simulations for 1000 degree preserving randomisations

1. HIvalidation.RData
Contains the:
	1.0. Transitivity and Averege Path Length (APL) simulations for 1000 degree preserving randomisations
	1.1. Graph and network HI objects
	1.2. The stratified random sample from the HI2 to validate the APL and transitivity in it.
