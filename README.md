decRiPPter (Data-driven Exploratory Class-independent RiPP TrackER) 

Alexander M. Kloosterman, Peter Cimermancic, Somayah S. Elsayed, Chao Du, Michalis Hadjithomas, Mohamed S. Donia, Michael A. Fischbach, Gilles P. van Wezel, and Marnix H. Medema.

Reference:
Integration of machine learning and pan-genomics expands the biosynthetic landscape of RiPP natural products
https://biorxiv.org/cgi/content/short/2020.05.19.104752v1

Description:

decRiPPter is a genome mining tool for detection of novel biosynthetic gene clusters (BGCs) of ribosomally synthesized and post-translationally modified peptides (RiPPs).
The tool prioritizes novelty at the cost of accuracy. As such, many of the BGCs that will come out of these results will not actually be RiPP BGCs. 
Rather, it functions as a platform for exploration and cluster prioritization, to detect RiPP BGCs that fall outside the scope of known RiPP classes. 

If you're more interested in accurate detection of RiPP BGCs of known classes, with a lower amount of novelty, 
we can refer you to some of the other excellent tools written for that purpose.
For example, see BAGEL4, RODEO, antiSMASH or RiPP-PRISM. 


Pipeline description

decRiPPter detects putative RiPP precursors using a Support Vector Machine (SVM) trained at the detection of RiPP precursors irrespective of RiPP subclass.
The genomic context of all candidate precursors is used to filter the results and narrow down to contain the features of interest. 
It then groups the remaining gene clusters together to form candidate families. Overlap with antiSMASH can also be analyzed, to remove known RiPP families.

decRiPPter is meant for the analysis of multiple closely related genomes simultaneously. 
From these genomes, it will infer information about the frequency of occurrence of genes (called to COG-score), to filter out household genes.
A large-scale analysis will also give a more accurate representation of the spread of a gene cluster family.

Description step by step:

1) SVM precursor detection
All encoded proteins with maximum length of 100 amino acid are analyzed with a pretrained SVM. In addition, all intergenic small open reading frames are also analyzed, even if these were not annotated.

2) COG analysis
The relative frequency of occurrence of each gene in all of the query genomes is determined next. Genes are grouped together if they have orthologue-like similarity to one another, the cutoff of which is based on the similarity of highly conserved genes in all genomes. See the publication for more details. 

3) Gene cluster formation
Operon-like gene clusters are formed around each candidate precursor. Only genes on the same strand as the precursor are used. 
Two methods are built in:
In the simple method, only intergenic distance is used as a cutoff
In the island method, genes are first fused within a small intergenic distance. The island may be further fused based on the difference of average COG scores.

4) Gene cluster annotation
Gene clusters are annotated with Pfam and TIGRFAM. 
A number of flanking genes is included to show additional context.
These domains are grouped into biosynthetic, transporter, regulator, peptidases and known RiPP categories.
The lists with these domains are found in the data/domains folder, and can be changed and expanded there.

5) Gene cluster filtering
Gene clusters are filtered based on passed filters
By default, outputs are generated for two different filter settings, called the mild filter and the strict filter.
See below and the publication for more details.

6) Gene cluster grouping
Gene clusters are grouped in two ways:
Precursor similarity (determined by blastp)
Jaccard index of protein domains found 

The resulting groups are further refined with the Markov Clustering Algorithm (MCL).
An additional grouping method is carried out when both precursors and protein domains are overlapping between two gene clusters. 
This last method is generally considered the most reliable (see publication). 


Protein fasta and genbank files will be generated as output for each operon. The formed gene cluster families are also written out as a network file that can be parsed by CytoScape.
In addition, for each of the filter settings, an index page will be generated which contains links to the families and graphical output of each formed gene cluster. 
On this index page, you can further filter the resulting gene clusters, based on gene cluster features. 


Installation guide:

decRiPPter is available as a commandline-tool for Linux.
As of now, it still runs in Python2 (Update to python3 is high on the todo list).

The easiest way to install it is to create a vritualenvironment

virtualenv decrippter -p python2

Then use pip to install these python packages. The following versions have been tested:

scikit-learn==0.11
biopython==1.76
scipy==1.2.3
matplotlib==2.2.5
networkx==2.2
numpy==1.16.6

In addition, make sure the following executables are in your $PATH variable.

blastp (from NCBI BLAST+, tested V2.6)
diamond (from DIAMOND, tested v0.9.31.132)
hmmsearch (from hmmer, tested V3.1b2)
mcl (from Markov Clustering Algorithm)

Optional:
prodigal (tested v2.6.3)

Later versions of these packages have not been tested, although they should work fine barring any changes in output files.
The only package with a version requirement is scikit-learn (v.011).

In addition, please download the latest versions of the Pfam and TIGRFAM databases,
from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/ and 
ftp://ftp.jcvi.org/pub/data/TIGRFAMs/, respectively.

Optional

Download and install antiSMASH V5 as specified on https://antismash.secondarymetabolites.org/#!/download.


Setting up the config file:
In the config file, let the variables pfam_db_path and tigrfam_db_path point to the Pfam and TIGRFAM databases, respectively.

Most settings in the config file can be given as arguments to the commandline application. 
Specifying them in the config file will save you the effort of typing them out every time.

Usage:

Quick usage:
Step 1
python genome_prep.py -o path/to/output -t taxid_to_download -i genomes_to_analyze PROJECT_NAME

Step 1.5 (Optional):
switch to antiSMASH environment
python antismash_wrapper.py -o path/to/output PROJECT_NAME
switch back to decRiPPter environment

Step 2
python genecluster_formation -o path/to/output PROJECT_NAME

Visual output can be found by opening the Index.html file in the output folder


Detailed usage:

Running the pipeline goes in two steps. 

Step 1. 

In the first step the genbank files are downloaded and/or parsed. Candidate precursors are detected and the COG scores are calculated.

usage: genome_prep.py [-h] [-o OUTPUTFOLDER] [-c CORES] [-i IN] [-t TAXID]
                      [-rg REUSE_GENOMES] [--run_prodigal {auto,never,always}]
                      [-p] [--load_pickles] [--store_cog] [--load_cog]
                      PROJECT_NAME

-c (CORES):             Number of processor cores to use.
-o (OUTPUTFOLDER):      The folder in which to create the project folder. The created project folder will be named after the project

1a) Input selection and genome downloading

- i (IN):               Point to a file or folder of files to be analyzed. Each file should correspond as one genome.
- t (TAXID):            Indicate the taxonomic identifier to download all genomes corresponding to that identifier. E.g., give "1883" to download all Streptomyces genomes from NCBI
--run_prodigal:         Reannotates downloaded files with prodigal. 
-rg (REUSE_GENOMES):    Used to reuse genome files already in the Genomes folder. Specify a file that contains a list of genomes to use exclusively, or use "all" to use all genomes in the Genome folder


Genomes in the folder specified by "-i" are analyzed in addition to the downloaded genomes (from "-t"). when redoing a run with the same genomes, simply use "-rg all" to reuse all genomes.
Either argument is optional, although at least one of the three arguments should be given.

By default all genomes under the given taxonomic identifier are downloaded. This can be further filtered in the config file.

Relevant config entries:

Which annotation to use when downloadig genomes from NCBI (genbank or refseq)

tax_file=refseq

When downloading, indicate which assembly levels are required, and which refseq_category they should belong (false for none)
Seperate multiple requirements with commas, without spaces. E.g. assembly_level_req=Scaffold,Complete Genome,Contig
Similar criteria can be passed to the required refseq category. E.g. representative genome,reference genome

assembly_level_req=false
refseq_category= false

1b) SVM

The cutoff is taken from the config file.

Relevant config entries:
SVM_cutoff=0.90

1c) COG calculation 

--store_cog:            Store the COG data seperately
--load_cog:             Load the COG data previously stored instead of calculating it

Useful when redoing runs

COG score calculations rely on a complete allVall BLAST, which is by default done with DIAMOND.
To prevent the results to be too large to be parsed, you can split up genomes randomly into subgroups.

Relevant config entries:
max_proteins=200000

The exact number ultimately relies on the amount of RAM of the system used to carry out the analysis.
For 8Gb of RAM, we estimate that 200,000 proteins is a realistic maximum.
The relationship is roughly quadratic: e.g. for ~200 Gb, you can go up to 1,000,000 proteins.

COG scores calculation relies on the identification of so-called trueCOGs: highly conserved genes present in all genomes of the subgroup analyzed.
If not enough trueCOGs are found, even smaller subgroups are formed. These groups must fulfill the criteria given here.

;Minimum genome group size
cog_min_group_size=5
;Minimum amount of truecogs
cog_min_truecogs_req=10

These groups are formed from the bottom up: starting with a random genome, genomes are randomly added. 
Each genome added will restrict which genes can be trueCOGs, since they must be conserved in all analyzed genomes.
Therefore, the number of trueCOGs will either stay the same or decrease with each added genome. 
The genomes can be added either in a random way (random), or each time the genome can be picked whose addition to the growing subgroup results in the smallest decrease in trueCOGs (best).
The latter option is more time-consuming than the former. 

cog_addmethod=random

After a subgroup is formed this way, a new seed genome is selected and the process is repeated. The resulting subgroups may overlap. 
The final goal is to let each genome be part of multiple subgroups that each fulfill the criteria, and cover as many other genomes as possible.
Even when new subgroups are created, they may not increase the number of genome pairs covered, as they were already covered by other subgroups. 
To prevent endless recreation of subgroups, the process is halted either after a certain amount of unsuccesfull iterations, or after a certain time (minutes). 

cog_stop_split_iterations=5
cog_stop_split_time=false

You can also set a minimum number of iterations, i.e. a minimum number of subgroups that would need to be formed at minimum before breaking off the group forming

cog_bottomsup_miniter = 3


Step 1.5)  Overlap with antiSMASH (optional)

Deactivate the current environment, and switch to the environment in which you can run antiSMASH.

usage: antismash_wrapper.py [-h] [-o OUTPUTFOLDER] [-c CORES] PROJECT NAME

Use the same arguments for "-o" and PROJECT_NAME as given to the genome_prep script. antiSMASH will be run in minimal mode on the genomes. 
After the script is done, switch back to the decRiPPter environment.

Make sure to set the parsing of antiSMASH files on in the config file. You can also specify the minimal overlap.
Results will be integrated in decRiPPter final output.

Relevant config entries:

antismash_parse=false
antismash_minimal_overlap=0.5
antismash_version=5

Step 2) Gene cluster formation, filtering and output generation

usage: Operon_formation.py [-h] [-o OUTPUTFOLDER] [-c CORES]
                           [--operon_method {simple,island}] [--skip-hmm]
                           [--load_operons]


Use the same arguments for "-o" and PROJECT_NAME as given to the genome_prep script.


2a) 

--operon_method:        Use from simple or island. If not given here, it will be read from the config file.

Further parameter can be given in the config file:

For simple gene cluster formation, use this as a distance cutoff.
simple_dist=750

For island gene cluster formation, set this as a maximum distance cutoff:
island_dist=750

For the formation of islands, the genes are fused using the following distance:
island_gene_dist=50

For fusing two islands, the difference in average COG scores is calculated. 
If the difference if below island_cog_cutoff + std_factor*(STD COG scores island 1 + STD COG scores island 2), the islands are fused. 

island_cog_cutoff=0.1
std_factor=1

2b) Annotation
Protein domain annotation with Pfam and TIGRFAM is carried out in the formed gene clusters plus flanking regions.


hmmsearch result files are kept, and are automatically reused. To overwrite the results, set the following config entry to always
run_hmm=auto

2c) Gene cluster filtering

Gene clusters are filtered based on the entries in the config file

gene_cluster_min_length=3
;Maximum average COG score of the gene cluster
gene_cluster_max_cog=0.25
;Minimum nr of genes encoding biosynthetic proteins, in the core gene cluster or anywhere (core or flanking region)
gene_cluster_domains_min_biosyn=2
gene_cluster_domains_min_biosyn_all=0

The filters can be specified to the core region or the entire gene cluster (including extension)
e.g. gene_cluster_domains_min_biosyn sets a requirement for the core region of the gene cluster, and gene_cluster_domains_min_biosyn_all for the core region plus the flanking genes

2c) Gene cluster grouping

Gene clusters are grouped in two ways, via precursor similarity or based on the Jaccard index of the detected protein domains.
MIBiG BGCs (1.4) are automatically added to the comparison

[clustering of BGCs]
;Cluster all BGCs found with any setting?
calc_networks=true
;Minimum group size
min_group_size=2

[precursor grouping]
;Whether or not to run the precursors allVall BLAST (will try to load previous results otherwise)
precursor_blast_run=true
;Cutoff for precursors to be grouped together
prec_min_ev=1
prec_min_bitscore=30
;Refine groups with mcl
precursor_mcl=true

[calculate_jaccard]
;Whether or not to group gene clusters with jaccard index of protein domains found in the core (will try to load previous results otherwise)
calculate_jaccard=true
;The cutoff to use
jaccard_cutoff=0.5
;Refine groups with mcl
jaccard_mcl=true

Finally, output is generated and can be found in the outputfolder

