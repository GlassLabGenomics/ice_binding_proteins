# Ice Binding Proteins (IBP) Project
This is a repository of scripts, documentation, and results for the computational pipeline of the ice binding protein (IBP) project: [**_Cold tolerance in the intertidal zone_**](https://github.com/GlassLabGenomics/ice_binding_proteins/wiki).

## Aims

We are developing a comprehensive screening method to reliably identify IBPs in intertidal invertebrates and gather information on their functional mechanisms, evolutionary history and ecological roles in Arctic and subarctic
regions. 

Our screening method is being developed as a combination of two approaches: 1) a screening pipeline algorithm, coupled with 2) ice-binding assays. This is a design that relies on interplay between computational and experimental techniques in an attempt to define and isolate a true ice-binding signal. The computational pipeline combs through current sequence and structure databases in search of proteins containing ice-binding motifs and their source organisms, and the experimental step verifies if ice-binding activity indeed occurs. We search for IBPs in two target species of sea stars, _Asterias amurensis_ (Northern Pacific seastar) and _Pisaster ochraceus_ (Ochre seastars), with the aim to expand the search to _Pycnopodia helianthoides_ (Sunflower seastar) and potentially other inverts. 

Our overall goal is to produce and test a method that can verify IBP existence and activity, and reliably predict novel IBPs. 

## Contents

1. Documentation: For a short description, see section on the screening pipeline below - for a more in-depth description see [wiki](https://github.com/GlassLabGenomics/ice_binding_proteins/wiki)
1. Code: See modules in [scripts](https://github.com/GlassLabGenomics/ice_binding_proteins/tree/main/scripts)
2. Results: See figures and descriptions in [results](https://github.com/GlassLabGenomics/ice_binding_proteins/results)

## Pipeline

There are two main questions we try to answer with our computational screening pipeline.  

**1. Which putative classes of IBPs are found in intertidal invertebrates (if any) and from which lineages do they arise?**

**2. What structures do these putative invertebrate IBPs adopt and how do they impact function?**

Our pipeline builds modules and connects existing tools to perform database querying with a verified set of confirmed IBP proteins (referred to as _seeds_ or _seed sequences_), searching through all available sequence and structural data. 
This querying is designed to be done from multiple angles to extract as much information as possible. Sequence querying thus occurs both at the gene and protein levels (i.e. searching against DNA/RNA databases and protein databases), and query results can analysed at whichever step of the taxonomic hierarchy, from large-scale clades to individual species. The emphasis with our computational pipeline is a usable way to filter through large amounts of data to identify putative IBPs and provide preliminary information regarding their evolution and prevalence in organismal groups to guide further targetted experimentation.

### Method

- preparation of input
    - prerequisite verification: target species data availability and quality (does a genome exist, are there annotations)
    - curating seed sequence dataset: adapting set from literature, what is already known to get hits
    - extraction of seed sequences: choice of database based on its characteristics such as sequence variety and known biases
- querying 
    - protaxovis tool
- output processing
    - taxonomic partitioning for viewing query results
    - individual sanity checks of hits vs target, alignment specs and overall reliability of query, done per protein
