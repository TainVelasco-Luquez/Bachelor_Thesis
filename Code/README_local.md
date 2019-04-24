# The files follows the same order as in the document:

1. `HIassembly.Rmd` Is the `R` markdown file, which contains the code used to download and cleaning up the interactions files for Protein-Protein Interactions (PPI), Mebolic-Signaling interactions (MetaboSignal) and regulatory interactions. This file also render the conversion to Entrez gene ID and, finally, the assembly of the brain-specific Human Interactome.
2. `HIassembly.html` is the html version of the HIassembly file, in which the code and substantial explanations for each step is present in an interactive file.
3. `HIvalidation.Rmd` Is the `R` markdown file, which contains the topological and functional validation procedures, used to show the biological relevance of the HI constructed. The .html file is still to come as, due to its size, its rendering is quite expensive, compuationally.
4. `RBDvalidation.Rmd` Is an `R` markdown file, which contains the code used to predict and validate, functionally and topologically, the putative REM sleep behavior disorder (RBD) module, according to the disease module hypothesis. The interactive .html file is still to come as, due to its size, its rendering is quite expensive, compuationally.
5. `Interactome_edge_proportion_per_source_comparison.R` Is an `R` script containing a comparison of four distinct interactome assemblies varying the the type of metabosignal edges (with or without compounds), regulatory edges (from Marbach et al. (2016) or from RegNet) and for PPI (APID level 2 or level 3). It is useful to decide which is the best combination of edge sources to gather into the final human interactome.
6. `RBD_robustness_APID_level2.R`, `RBD_robustness_APID_level3.R`, `RBD_robustness_brain_adult_2.5.R`, `RBD_robustness_brain_adult_greater_mean.R`, `RBD_robustness_metabosignal_ready.R`,`RBD_robustness_RegNetwork.R` are all `R` scripts finding the RBD disease module with diferent input networks (described in the file name).
7. `Regulatory_exploration.R` is an `R` script exploring some basic features of the Marbach et al (2016) brain regulatory information.
8. `Robustness_powerlaw_brain_adult.R` is an `R` script assessing the robustness of the powerlaw exponent when deleting a given percentage of the interactions in the network.
8. `DIAMOnD_stop_criterion_comparison_among_diferent_datasets.R`is an `R` script comparing the DIAMOnD stop criterion among all sources of information and specially between Marbach et al (2016) and RegNet.
9. `interactome_with_brain_adult_2.5.R` is an `R` script which finds the RBD disease module in the Human Interactome with only the 2.5% of regulatory interactions from Marbach et al (2016). I.e. the HI contains 2.5% Marbach + Full metabosignal + Full APID level 3.

## Notes
If the absolute path to a file inside the abovementioned scripts is `"~/Documents/Thessis/filename"`, then the given file must be in the relative path `"Data/Individual/filename"`.


## Todo list
* Parallelise the assembly of the `EG2GOBP` object.
