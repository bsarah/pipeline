![Pipeline](http://www.bioinf.uni-leipzig.de/~bsarah/tRNAinsertion.png "")

##Pipeline

This is the README of our pipeline (still need to find a name).
Use this is a tool to detect homologous genetic elements in different species
based on a model or a list of the genetic elements and a multiple sequence
alignment of the target species. You can choose if secondary structure
informations should be included in the analysis or not.

#Program use:

perl main.pl -o OutputFolder -g GenomesFolder -m MafFolder -f refspecies
[-s seqSimilarity] [-t strucSimilarity] [-q threshold] {-c CM | -l genelist}
-y python -e perl -r R -i infernal

[] means optional

{} means one of the options should be given

all other parameters are required

#Parameter explanations:

-o Outputfolder:
   please specify where the Output files should be created.

-g GenomesFolder:
   path to folder where genomes are located. The names of the genome files must
   be the same as the species' name in the MSA. The genome files should have
   fasta format and can be gzipped.

-m MafFolder:
   path to the folder where multiple sequence alignments in maf format are
   located. The names of the species in the maf files have to be the same as
   the names of the genome files. Format should be: .maf.gz, .maf.Z, .maf.bz2

-f refspecies:
   the maf files were created by using one species as the reference species.
   Please provide the name of the reference species as it can be found in the
   multiple sequence alignment and as a genome file

-s seqSimilarity and -t strucSimilarity:
   when testing for homology, only sequences with a sequence similarity above
   seqSimilarity and a secondary structure similarity above strucSimilarity
   will be considered as homologous whereas the value is measured as percentage
   , thus 1 means the sequences are the same. Default for both : 0.9; if for
   any of the
   two parameter -1 is given, this information will be skipped anthe analysis
   will only be based on one of both, sequence or structure information.

-q threshold:
   remove a percentage of the lowest scoring blocks based on the MAF scores,
   whereas the value is between 0 and 100

-c CM or -l genelist:
   either give a covariance model (cm file) for the gene of interest or give a
   folder containing a list of genes for every species (format: .bed, .bed.gz,
   .bed.Z, .bed.bz2).
   The gene list should contain one line per gene with tab separated elements
   as follows:

   chromosome     species     startCoord     endCoord     '+'or'-'    sequence
   secondary_structure(optional).

   If either secondary structure or sequence are not given or should not be
   taken
   into account, please use the parameters -s and -t and set the corresponding
   one to -1.

-w for inclusion E-value threshold or -b for an inclusion bit score threshold
   as parameter for the infernal run:

   As described in the infernal manual (to be found here: eddylab.org/infernal)
   , inclusion thresholds control which hits are considered to be significant.
   In case you want to include more putative sequences into our pipeline as
   significant hits based on the given
   covariance model, please either increase the e-value threshold using
   parameter -w (default = 0.01) or decrease
   the bitscore threshold with parameter -b (this value is usually not set
   when running infernal with default
   parameters). Infernal will consider sequences significant if their
   e-value <= the e-value threshold or their
   bit score >= the bit score threshold. The README_output file will explain
   where to find the infernal output in
   order to see which sequences were included in the analysis.

-y python:
   please specify the path to where python is installed. Version should
   be >= 3.0

-e perl:
   please specify the path to where perl is installed. Version should be >= 5.0

-r R:
   please specify the path to where R is installed. Version should be >= 3.2.0

-i Infernal:
   please specify the path to where Infernal is installed. Version should
   be >= 1.1.1


#Further parameters:

choose -k for contact information and citation, -v to see the
version information and -h to display this help message.
