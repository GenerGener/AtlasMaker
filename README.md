# AtlasMaker

The goal of this repo will be to abstract methods from hiv-atlas-maker.py development to create species-agnostic method for reference-guided mRNA modeling.

Note that because of U3-R ideosyncracies and how users might supply their input and reference genomes, .... The latest implentation works best if the user supplies a genome that is longer than the mRNA's which were used to inform HIV-1 M NL4-3 mRNA (HIV-1 Atlas V1.0). The unspliced HIV-1 M NL4-3 mRNA model [MZ242719.1](https://www.ncbi.nlm.nih.gov/nuccore/MZ242719.1/) is 9173 bases long. Other types of loci which might be amenable to the approache here include repetitive elements, whether endogenous or exogenous (e.g., retroviruses), which leverage direct long terminal repeats in their propagation. 

Related repo: [AtlasMatcher](https://github.com/GenerGener/AtlasMatcher)

