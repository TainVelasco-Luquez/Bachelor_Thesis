# Thesis
R markdown files and data sets required for make reproducible my graduate thesis entitled "On the REM sleep behavior disorder module: A network medicine approach", readable in https://www.overleaf.com/read/gktqggmgnjzp.

The files follows the same order as in the thesis:

    0. HIassembly.Rmd Is the R markdown file, which contains the code used to download and cleaning up the interactions files for Protein-Protein Interactions (PPI), Mebolic-Signaling interactions (MetaboSignal) and regulatory interactions. This file also render the conversion to Entrez gene ID and, finally, the assembly of the brain-specific Human Interactome.
    
    1. HIassembly.html is the html version of the HIassembly file, in which the code and substantial explanations for each step is present in an interactive file.
    
    2. HIvalidation.Rmd Is the R markdown file, which contains the topological and functional validation procedures, used to show the biological relevance of the HI constructed. The .html file is still to come as, due to its size, its rendering is quite expensive, compuationally.
    
    3. RBDvalidation.Rmd Is an R markdown file, which contains the code used to predict and validate, functionally and topologically, the putative REM sleep behavior disorder (RBD) module, according to the disease module hypothesis. The interactive .html file is still to come as, due to its size, its rendering is quite expensive, compuationally.
    
I hope readers find both the document and the code enjoyable. Contributions and critics are more than welcome from anybody interested, feel free to contact me if you find somethig wrong. 

To do:
1. Replace the absolute paths for relative ones. This wiil make the downloaded directory automatically selfcontained and easy accesible from the scripts
2. Optimise the code
