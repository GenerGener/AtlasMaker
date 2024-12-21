# AtlasMaker

The goal of this repo will be to abstract methods from hiv-atlas-maker.py development to create species-agnostic method for reference-guided mRNA modeling.

Note that because of U3-R ideosyncracies and how users might supply their input and reference genomes, .... The latest implentation works best if the user supplies a genome that is longer than the mRNA's which were used to inform HIV-1 M NL4-3 mRNA (HIV-1 Atlas V1.0). The unspliced HIV-1 M NL4-3 mRNA model [MZ242719.1](https://www.ncbi.nlm.nih.gov/nuccore/MZ242719.1/) is 9173 bases long. Other types of loci which might be amenable to the approache here include repetitive elements, whether endogenous or exogenous (e.g., retroviruses), which leverage direct long terminal repeats in their propagation. 

Related repo: [AtlasMatcher](https://github.com/GenerGener/AtlasMatcher)

# Notes from 20241219

Research question: *Does isoform-informed alignment increase sequencing depth, reference coverage?*

Standard approach maps to species-specific ref. 

While viruses are sometimes geographically associated, because of globalization, it isn’t safe to assume a given person’s HIV’s subtype. HIV is often assumed to be HIV-1.

Current bioinformatics pipelines often attempt to align sequencing reads to either HIV-1-only or human+HIV-1. A single HIV-1 “species” reference is often used. (Cite my HIV plasmid preprint). This risks missing out on information at the per-base sequence level if the HIV in the sample is divergent from the reference used. 

If multiple references are used, a fine balancing act between matching reference sequence and bypassing intrasample genome heterogeneity ensues. Depending on the aligner’s settings, and how HIV sequences are provided to the aligner, reads might be split across HIV “chromosomes”(sometimes called haplotypes, or whole-genome genotypes). Add to this picture the existence of circulating recombinant forms of viruses, and it becomes tricky to try to pick a more specific reference genome beforehand.

If time permits, a sample can be compared to a set of representative HIV genomes, and a best match determined by selecting the one with least divergence, answering the question “what is this sample’s HIV’s closest match?”

If this is done, this might “help” the aligner then recover the most of that HIV’s reads from a given sequencing experiment.

However, there exists another layer of HIV-seq heterogeneity that is independent of per-base information: HIV mRNA splicing. Splicing effectively breaks up sequence. If the aligner is not splice aware, splice reads will remain unaligned/unmapped. 

## Figure with worse to least mapping

Natural HIV reads aligned to species reference
Natural HIV reads compared to a set of HIV references. Best match used for sample-specific mapping.
Natural HIV reads aligned to species reference transcript models
Natural HIV reads aligned to sample-specific imputed transcript models (proposed)
Natural HIV reads aligned to sample-specific transcript models (impractical).

Note the natural HIV’s genome can be assembled de novo if both seq depth and reference coverage are sufficient. Reads can be enriched from sample background by selecting reads that exhibit homology to HIV reference genomes per the spectrum above. Iteration after each step can reassure the investigator about match quality.



