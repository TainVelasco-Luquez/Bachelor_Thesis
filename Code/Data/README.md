Description of the `R` objects inside each file linked to the `HIassembly.Rmd`, `HIvalidation.Rmd` and `RBDmodule.Rmd`, in order to make the rendering of such `.Rmd` faster:

1. `RBDmodule.RData` contains the most computationally expensive objects created in the file `RBDmodule.Rmd`:
    	1. igraph and network RBDmodule objects
	2. ERGM fit
	3. NEAT
	4. Transitivity and Averege Path Length (APL) simulations for 1000 degree preserving randomisations
2. `HIvalidation.RData` contains:
	1. Transitivity and Averege Path Length (APL) simulations for 1000 degree preserving randomisations
	2. Graph and network HI objects
	3. The stratified random sample from the HI2 to validate the APL and transitivity in it.
