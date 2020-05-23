# decRiPPter (Data-driven Exploratory Class-independent RiPP TrackER) 

Alexander M. Kloosterman, Peter Cimermancic, Somayah S. Elsayed, Chao Du, Michalis Hadjithomas, Mohamed S. Donia, Michael A. Fischbach, Gilles P. van Wezel, and Marnix H. Medema.

Reference:
Integration of machine learning and pan-genomics expands the biosynthetic landscape of RiPP natural products
Preprint: https://biorxiv.org/cgi/content/short/2020.05.19.104752v1

Sample output:
Large scale Streptomyces analysis (1,295 genomes)
[Mild filter] (http://www.bioinformatics.nl/~medem005/decRiPPter_mild/index.html)
[Strict filter] (http://www.bioinformatics.nl/~medem005/decRiPPter_strict/index.html)

## Description:

decRiPPter is a genome mining tool for detection of novel biosynthetic gene clusters (BGCs) of ribosomally synthesized and post-translationally modified peptides (RiPPs).

decRiPPter functions as a platform for the exploration and prioritization of candidate RiPP BGCs. **It prioritizes novelty at the cost of accuracy**. As such, many of the BGCs that will come out of these results may not actually be RiPP BGCs. However, **the resuling BGCs are not restricted to specific genetic markers, and may be novel RiPP BGCs as a result**. To help interpret the results, it allows for extensive user-defined filtering of the candidate RiPP BGCs, to detect RiPP BGCs that fall outside the scope of known RiPP classes. 

If you're more interested in highly accurate detection of RiPP BGCs, staying within the bounds of known RiPP subclasses, 
there are some excellent tools written for that purpose.
For example, see [BAGEL4](http://bagel4.molgenrug.nl/), [RODEO](http://ripp.rodeo/), [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) or [RiPP-PRISM](https://github.com/magarveylab/prism-releases). 


decRiPPter detects putative RiPP precursors using a Support Vector Machine (SVM) trained for the detection of RiPP precursors irrespective of RiPP subclass.
The genomic context of all candidate precursors is used to filter the results and narrow down to contain the features of interest. 
It then groups the remaining gene clusters together to form candidate families. Overlap with antiSMASH can also be analyzed, to remove known RiPP families.


**decRiPPter is meant for the analysis of multiple closely related genomes simultaneously.**
From these genomes, it will infer information about the frequency of occurrence of genes (called to COG-score), to filter out household genes.
Analyzing more genomes simultaneously allows gives a representation of the spread of a gene cluster family, which can help in the filtering process.


### Pipeline description

#### 1) SVM precursor detection
All encoded proteins with maximum length of 100 amino acid are analyzed with a pretrained SVM. In addition, all intergenic small open reading frames are also analyzed, even if these were not annotated.


#### 2) COG analysis
The relative frequency of occurrence of each gene in all of the query genomes is determined next. Genes are grouped together if they have orthologue-like similarity to one another, the cutoff of is based on the similarity of highly conserved genes in all genomes. See the publication for more details.


#### 3) Gene cluster formation
Operon-like gene clusters are formed around each candidate precursor. Only genes on the same strand as the precursor are used. 

Two methods are built in:
In the simple method, only intergenic distance is used as a cutoff
In the island method, genes are first fused within a small intergenic distance. The island may be further fused based on the difference of average COG scores


#### 4) Gene cluster annotation
Gene clusters are annotated with Pfam and TIGRFAM. 
A number of flanking genes is included to show additional context.
These domains are grouped into biosynthetic, transporter, regulator, peptidases and known RiPP categories.
The lists with these domains are found in the data/domains folder, and can be changed and expanded there.


#### 5) Gene cluster filtering
Gene clusters are filtered based on passed filters
By default, outputs are generated for two different filter settings, called the mild filter and the strict filter.
See below and the publication for more details.


#### 6) Gene cluster grouping
Gene clusters are grouped in two ways:
Precursor similarity (determined by blastp)
Jaccard index of protein domains found 


The resulting groups are further refined with the Markov Clustering Algorithm (MCL).
An additional grouping method is carried out when both precursors and protein domains are overlapping between two gene clusters. 
This last method is considered the most reliable by the authors (see preprint). 

#### 7) Output generation
Protein fasta and genbank files will be generated as output for each operon. The formed gene cluster families are also written out as a network file that can be parsed by CytoScape.
In addition, for each of the filter settings, an index page will be generated which contains links to the families and graphical output of each formed gene cluster. 
On this index page, you can further filter the resulting gene clusters, based on gene cluster features. 


## Installation guide:

decRiPPter is available as a commandline-tool for Linux.
As of now, it still runs in Python2 (Update to python3 is high on the todo list).

#### 1) Clone the environment locally

```git clone ```

#### 2) Setup the environment

The easiest way to install it is to create a vritual environment

```virtualenv decrippter -p $(which python2) decrippter```

Then use pip to install these python packages. The following versions have been tested:

scikit-learn==0.11

biopython==1.76

scipy==1.2.3

matplotlib==2.2.5

networkx==2.2

numpy==1.16.6

These can be installed using the ```requirements.txt``` file

```pip install -r requirements.txt```

In addition, make sure the following executables are in your ```$PATH``` variable.

```blastp``` (from [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK52640/), tested V2.6)

```diamond``` (from [DIAMOND](https://github.com/bbuchfink/diamond), tested v0.9.31.132)

```hmmsearch``` (from [hmmer](http://hmmer.org/), tested V3.1b2)

```mcl``` (from the [Markov Clustering Algorithm](https://micans.org/mcl/))



In addition, please download the latest versions of the [Pfam](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/) and [TIGRFAM](ftp://ftp.jcvi.org/pub/data/TIGRFAMs/) databases,


Later versions of these packages have not been tested, although they should work fine barring any changes in output files.
The only package with a version requirement is scikit-learn (v.011).


**Optional**

```prodigal``` (from [prodigal](https://github.com/hyattpd/Prodigal), tested v2.6.3)

Genome (re)annotation is built in with decRiPPter via prodigal, although it is not a strict requirement.


```antismash``` (from [antiSMASH V5](https://antismash.secondarymetabolites.org/#!/download))

Download and install antiSMASH V5 as specified in it's own environment. 

#### 3) Setup the config file:
In the config file, let the variables ```pfam_db_path``` and ```tigrfam_db_path``` point to the Pfam and TIGRFAM databases, respectively.


## Usage:


### Quick start:


#### Step 1

    python genome_prep.py -o path/to/output -t taxid_to_download -i genomes_to_analyze PROJECT NAME


#### Step 1.5 (Optional):

Deactivate the decRiPPter environment

```deactivate```

Switch to antiSMASH environment, e.g. if antiSMASH is installed via conda:

```conda activate antismash``` 

Run the antiSMASH wrapper; use the same arguments for ```-o``` and ```PROJECT NAME```

    python antismash_wrapper.py -o path/to/output PROJECT NAME

Switch back to decRiPPter environment

    source /path/to/decrippter_env/bin/activate


#### Step 2

    python genecluster_formation -o path/to/output PROJECT_NAME

Visual output can be found by opening the Index.html file in the output folder


### Detailed usage:

Running the pipeline goes in two (optionally three) steps. 

#### Step 1. 

In the first step the genbank files are downloaded and/or parsed. Candidate precursors are detected and the COG scores are calculated.

```usage: genome_prep.py [-h] [-o OUTPUTFOLDER] [-c CORES] [-i IN] [-t TAXID]```
                      ```[-rg REUSE_GENOMES] [--run_prodigal {auto,never,always}]```
                      ```[-p] [--load_pickles] [--store_cog] [--load_cog]```
                      ```PROJECT NAME```

```-c``` (CORES):             Number of processor cores to use.

```-o``` (OUTPUTFOLDER):      The folder in which to create the project folder. The created project folder will be named after the project

##### 1a) Input selection and genome downloading

```-i``` (IN):               Point to a file or folder of files to be analyzed. Each file should correspond to one genome.

```-t``` (TAXID):            Indicate the taxonomic identifier to download all genomes corresponding to that identifier. E.g., give "1883" to download all Streptomyces genomes from NCBI. By default all genomes under the given taxonomic identifier are downloaded. Additional requirements for downloading genomes can be set in the config file. 

```--run_prodigal (PRODIGAL_MODE)```:         Reannotates downloaded files with prodigal. When set to ```never```, the program will only download genbank files and parse these. When set to ```always```, the program will only download DNA fasta files, and annotate these to create a very basic genbank file. When set to ```auto``` (default), the pipeline will download/use genbank files when available. If these are not found, it will download DNA fasta files instead and annotate them.

```-rg``` (REUSE_GENOMES):    Used to reuse genome files already in the Genomes folder. Specify a file that contains a list of genomes to use exclusively, or use "all" to use all genomes in the Genome folder


Genomes in the folder specified by ```-i``` are analyzed in addition to the downloaded genomes (from ```-t```). when redoing a run in this step with the same genomes, simply use ```-rg all``` to reuse all genomes. Either argument is optional, although at least one of the three arguments should be given.

Additional details for downloading can be set in the config file.

Which annotation to use when downloadig genomes from NCBI (genbank or refseq)

```tax_file=refseq```

When downloading, indicate which assembly levels are required, and which refseq_category they should belong (false for none)
Seperate multiple requirements with commas, without spaces. E.g. assembly_level_req=Scaffold,Complete Genome,Contig
Similar criteria can be passed to the required refseq category. E.g. representative genome,reference genome

```assembly_level_req=false```

```refseq_category= false```

##### 1b) SVM analysis of precursors


The cutoff can be set in config file (0-1). 

Relevant config entries:

```SVM_cutoff=0.90```

You can additionally set the maximum overlap between annotated genes when detecting small open reading frames in intergenic regions.

```maximum_overlap_smorf=20```

##### 1c) COG calculation 

```--store_cog```:            Store the COG data seperately

```--load_cog```:             Load the COG data previously stored instead of calculating it

Useful when redoing runs, as it foregoes the entire COG calculation. Can be used when rerunning with a subset of the original genomes. 

COG score calculations rely on a complete allVall BLAST, which is by default run by DIAMOND.
To prevent the results to be too large to be parsed, you can set an upper limit in the config file for how many proteins are allowed in the groups of genomes analyzed. If the number of proteins exceeds this number, the genomes will be randomly and evenly split until this is no longer the case.

```max_proteins=200000```

The exact number relies on the amount of RAM of the system used to carry out the analysis.
For 8Gb of RAM, 200,000 proteins is a realistic maximum in our experience.
The relationship scales roughly quadratic: e.g. for ~200 Gb, you can go up to 1,000,000 proteins.

COG scores calculation relies on the identification of so-called trueCOGs: highly conserved genes present in all genomes of the subgroup analyzed.
If not enough trueCOGs are found, even smaller subgroups are formed. These groups must fulfill the criteria given here.

```cog_min_group_size=5```

```cog_min_truecogs_req=10```

Always make sure that you're analyzing enough genomes to fulfill the minimum group requirement.

The final goal of forming all possible subgroups is to let each genome be part of a collection subgroups that each fulfill the criteria set above. These subgroups are allowed to be overlapping, so that ideally, each genome pair can be formed via any of the groups formed.

However, we can not possibly form all possible subgroups of genomes and test which are viable, since this process is incredibly time-consuming. Therfore, the following method is used.

These groups are formed from the bottom up: starting with a random genome, genomes are added one at a time as long as the groups fulfill the criteria.  
Each genome added will restrict which genes can be trueCOGs, since they must be highly conserved in the newly added genome too. Therefore, the number of trueCOGs will either stay the same or decrease with each added genome, until the number of trueCOGs reaches the criteria set or all genomes have been joined.

The genomes can be added either in a random way (random), or each time the genome can be picked whose addition to the growing subgroup results in the smallest decrease in trueCOGs (best).
The latter option is more time-consuming than the former, as for each addition it needs to be calculated which genome addition would result in the smallest decrease. However, it may result in a higher number of trueCOGs, especially in small to intermediate-sized groups.

Choose random or best

```cog_addmethod=random```

After group formation, a new genome is selected to start a new group. This genome is selected from the genomes that have the lowest number of genome pairs set collectively in all their groups. 

As long as the addition of new subgroups increases the quality of the subgroups (i.e. more genome pairs are formed), the process continues. The process is halted if the genome pairs have no longer improved (e.g. no new genome pairs were added) after a certain amount of iterations, or after a certain time (minutes). 

```cog_stop_split_iterations=5```

```cog_stop_split_time=false```

When a genome pair is covered by two subgroups, the larger of the two will always be chosen.

In addition, the process stops if all genomes are at least part of one subgroup and the minimum number of iterations has passed.

```cog_bottomup_min_iter=10```


#### Step 1.5)  Overlap with antiSMASH (optional)

Deactivate the current environment, and switch to the environment in which you can run antiSMASH (see quick usage). 

usage: antismash_wrapper.py [-h] [-o OUTPUTFOLDER] [-c CORES] PROJECT NAME

Use the same arguments for ```-o``` and ```PROJECT NAME``` as given to the genome_prep script. antiSMASH will be run in minimal mode on the genomes. The result for each genome can be found in the genome folder. 
After the script is done, switch back to the decRiPPter environment.

When doing this step, make sure to set the parsing of antiSMASH files on in the config file. Results will be integrated in decRiPPter final output, and can be used to filter (see below).

```antismash_parse=false```

You can also specify the minimal fraction a decRiPPter gene cluster must be covered by an antiSMASH gene cluster to be considered an overlap. 

```antismash_minimal_overlap=0.5```

#### Step 2) Gene cluster formation, filtering and output generation

```usage: gene_cluster_formation.py [-h] [-o OUTPUTFOLDER] [-c CORES]```
                             ```[--operon_method {simple,island}] [--skip-hmm]```
                             ```[--load_operons] PROJECT NAME```


Use the same arguments for ```-o``` and ```PROJECT_NAME``` as given to the genome_prep script.

If you want to rerun the gene cluster formation with different settings, but want to keep your old results, rename the folder ```Output``` to something else, and rerun the script. 

2a) gene cluster formation

```--operon_method```:        Use from simple or island. If not given here, it will be read from the config file.

Further parameters can be given in the config file.

For simple gene cluster formation, use this as a distance cutoff.

```simple_dist=750```

For island gene cluster formation, the genes are fused into islands using the following distance:

```island_gene_dist=50```

For island gene cluster formation, set this as a maximum distance cutoff between islands

```island_dist=750```

For fusing two islands, the difference in average COG scores is calculated. 
If the difference if below ```island_cog_cutoff + std_factor*(STD COG scores island 1 + STD COG scores island 2```, the islands are fused. 

    island_cog_cutoff=0.1
    std_factor=1

##### 2b) Extension and annotation
The core, operon-like region of the gene cluster is extended with a number of flanking genes, which are added regardless of strand or distance. These are mostly used for downstream filtering and visualization. You can set the number to be added on each side in the config file.

```genecluster_extension=5```

Protein domain annotation with Pfam and TIGRFAM is carried out in the formed gene clusters plus flanking regions. Make sure the paths are properly set.

    pfam_db_path=/path/to/Pfam-A.hmm
    tigrfam_db_path=/path/to/TIGRFAMs_15.0_HMM.LIB

hmmsearch result files are kept, and are automatically reused in runs when found. To overwrite the results, set the following config entry to always. 

```run_hmm=auto```

Based on the domains found, genes can be annotated as encoding for a biosynthetic, transporter, regulator, peptidase or known RiPP (kripp) product. To see which domains are used for this annotation, check out the ```data/domains``` folder, downloaded with decRiPPter.

You can leave out any domain in the list by adding a # sign for it, or add new Pfam/TIGRFAM domains to your liking.

##### 2c) Gene cluster filtering

Gene clusters are filtered based on the entries in the config file

Minimum number of genes in the core region of the gene cluster

```gene_cluster_min_length=3```

Maximum average COG score of the gene cluster (core region only)

```gene_cluster_max_cog=0.25```

When using the optional antiSMASH annotation, you can filter these out here.

```filter_antismash=false```

Additional requirements are set for a minimum number of genes of each of the categories. 
E.g. Minimum number of genes encoding biosynthetic proteins, in the core region or in the core with flanks .

    gene_cluster_domains_min_biosyn=2
    gene_cluster_domains_min_biosyn_all=0

The filters can be specified to the core region or the entire gene cluster (including extension).
E.g. ```gene_cluster_domains_min_biosyn``` sets a requirement for the core region of the gene cluster, and ```gene_cluster_domains_min_biosyn_all``` for the core region plus the flanking genes



##### 2c) Gene cluster grouping

Gene clusters are grouped in two ways, via precursor similarity or based on the Jaccard index of the detected protein domains. For this purpose, only the core regions of the gene clusters are used. 
MIBiG RiPP BGCs (V1.4) are downloaded in a preformatted file and automatically added to the comparison.

###### Clustering of BGCs
Minimum group size

```min_group_size=2```

###### Precursor grouping
Cutoff for precursors to be grouped together, in E-value and bitscore, respectively.

    prec_min_ev=1
    prec_min_bitscore=30

Refine groups with mcl; weights are the percentage ID similarities

```precursor_mcl=true```

###### Domain-based grouping

The cutoff to use

```jaccard_cutoff=0.5```

Refine groups with mcl; weights are the Jaccard indeces

```jaccard_mcl=true```

Finally, output is generated and can be found in the folder ```Output```.

