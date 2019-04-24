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
<<<<<<< HEAD
3. `RBDmodule_NEAToutput.RData` contains:
  1. `RBDrich` is the NEAT (Network Enrichment Analysis Test) output when enriching RBD protomodule genes with the GO terms selected as described in `RBDmodule.Rmd`.
  2. `validationGOBP` is the NEAT output when enriching RBD protomodule genes plus the DIAMOND selected nodes per iteration, as described in `RBDmodule.Rmd`.
=======
3. `diamond_validation_neat.RData` contains:
    1. EG2GOBP is a data frame mapping GO IDs to EG IDs of the genes present in the RBD proto module
    2. RBDrich is a data frame from NEAT output with the network enrochment for the RBD proto module
    3. RBDrich10 is a data frame with 10 rows corresponding to the 10 most significantly enriched GO IDs in the RBD proto module
    4. validationGOBP100 a data frame with the output from the loop described in the file `Code/RBDmodule.Rmd` for the first 100 iterations of DIAMOnD
    5. validationGOBP500 a data frame with the output from the loop described in the file `Code/RBDmodule.Rmd` for the 500 iterations of DIAMOnD
>>>>>>> 68218bd5c97078d3282d429ebdef9feb711fbe0c
4. `export_Mon_Mar_18_13_56_38_UTC_2019.csv` is the RegNetwork downloaded files.
5. `HIassembly.RData` contains three objects: metabosignal, ppi and regulatory all ready interactions imported from the files with the same names in the Individual folder
6. `Data/interactome.model.ergm.RData` contains fourth files: 1. `interactomeGraph.model` from chunk ergm_1 in `HI_validation.Rmd`, 2. `interactomeGraph.model.gof.deg` from chunk goodness_of_fit_degree in `HI_validation.Rmd`, 3. `interactomeGraph.model.gof.edge` from chunk goodness_of_fit_edges in `HI_validation.Rmd`, 4. `interactomeGraph.model.gof.dis` from chunk goodness_of_fit_geodesic in `HI_validation.Rmd`.
