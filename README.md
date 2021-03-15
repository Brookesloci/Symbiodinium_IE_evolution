# Dinoflagellate Intron Evolution

This repository harbors data and documented scripts corresponding to the analyses conducted in the following manuscript:


**Scott William Roy, Landen Gozashti, Bradley A. Bowser and Brooke N. Weinstein, Graham E. Larue (2020). Massive intron gain in the most intron-rich eukaryotes is driven by introner transposable elements of unprecedented diversity and flexibility. BioRxiv https://doi.org/10.1101/2020.10.14.339549**


### Ancestral intron reconstruction

The folder, `intron_models`, contains the alignment files and phylogenetic tree file to reproduce the ancestral intron reconstruction in the paper using the graphical user interface in [Malin](http://www.iro.umontreal.ca/~csuros/introns/malin/) along with results from the analyses.


### Identifying and analyzing introner families

The folder, `ILEs`, contains fasta files of sequences for each ILE family initially identified in each considered species, along with consensus sequences for each ILE family and candidate ILE sequences found in intergenic regions. The folder also contains the custom perl scripts used to initially identify ILEs from annotated genomes and identify phase bias with regard to new insertions. Furthermore, we provide as series of python scripts used to construct consensus sequences for ILE families and retrieve candidate insertions in intergenic regions. 

#### ILE data for species

Files ccontaining initially identified ILEs, ILE candidates from intergenic regions, and ILE consensus sequences for each species are contained within a corresponding tar ball. For example, Symbiodinium A data is contained within `ILEs/symA.tar.gz`. 
