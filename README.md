# Ice Binding Proteins (IBP) Project
This is a repository of scripts, documentation, and results for the computational pipeline of the ice binding protein (IBP) project: [**_Cold tolerance in the intertidal zone_**](https://github.com/GlassLabGenomics/ice_binding_proteins/wiki).

## Aims

We are developing a comprehensive screening method to reliably identify IBPs in intertidal invertebrates and gather information on their functional mechanisms, evolutionary history and ecological roles in Arctic and subarctic
regions. 

Our screening method is being developed as a combination of two approaches: 1) a screening pipeline algorithm, coupled with 2) ice-binding assays. This is a design that relies on interplay between computational and experimental techniques in an attempt to define and isolate a true ice-binding signal. The computational pipeline combs through current sequence and structure databases in search of proteins containing ice-binding motifs and their source organisms, and the experimental step verifies if ice-binding activity indeed occurs. We search for IBPs in two target species of sea stars, _Asterias amurensis_ (Northern Pacific seastar) and _Pisaster ochraceus_ (Ochre seastars), with the aim to expand to _Pycnopodia helianthoides_ (Sunflower seastar) and potentially other inverts. 

Our overall goal is to produce and test a method that can verify IBP existence and activity, and reliably predict novel IBPs. 

## Contents

1. Documentation: For a short description, see section on the screening pipeline below - for an expanded description see [wiki](https://github.com/GlassLabGenomics/ice_binding_proteins/wiki)
1. Code: See modules in [scripts](https://github.com/GlassLabGenomics/ice_binding_proteins/tree/main/scripts)
2. Results: See figures and descriptions in [results](https://github.com/GlassLabGenomics/ice_binding_proteins/results)

## Pipeline

There are two main questions we try to answer with our computational screening pipeline. We do this through database querying with a verified set of confirmed IBP proteins (referred to as _seeds_), searching through all available sequence and structural data. 

**1. Which classes of IBPs are found in intertidal invertebrates and from which lineages do they arise?**

**2. What structures do the IBPs adopt and how do they impact function?**
