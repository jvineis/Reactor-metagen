# Reactor-metagen
## Documentation of the metagenomic analysis including MAG collections, MAG annotation, phylogenetics, phylogenomics, and pangenomics. The quality filtered reads used as the input for the assembly and all mapping steps are found here (figshare link in progress). The metadata for each of the samples is contained in the supplemental table (S1) in the manuscript (link to mscript).

1. Assemble the reads. You will need all of the metagenomic datasets and the file "quality-filtered-data-list.txt" in a single directory to run the following command. If you have a SLURM management system, please contact me for alternate methods to run this step through an array submission.
     
    #!/bin/bash
    module load megahit/1.0.6
    for i in `cat quality-filtered-data-list.txt`
    do
       gunzip "$i".gz
       megahit -r "$i" -o "$i"-MEGAHIT
       gzip "$i"
    done

2. Map the reads.
    #!/bin/bash
    for i in `cat ../jhv-samples.txt`
    do
        bowtie2 -x anvi-contigs -f -U ../"$i" -S "$i".sam
        samtools view -bS -F 4 "$i".sam > "$i"-RAW.bam
        rm -rf "$i".sam
        anvi-init-bam "$i"-RAW.bam -o "$i".bam
        rm -rf "$i"-RAW.bam
    done

3. Build contigs db.

    #!/usr/bin/bash
    anvi-script-reformat-fasta final.contigs.fa -o anvi-contigs.fa -l 2000 --         simplify-names
    bowtie2-build anvi-contigs.fa anvi-contigs
    anvi-gen-contigs-database -f anvi-contigs.fa -o	anvi-contigs.db
    anvi-run-hmms -c anvi-contigs.db
     
4. profile

    #!/bin/bash
    for i in `cat ../jhv-samples.txt`
    do
        clusterize anvi-profile -c anvi-contigs.db -T 5 -i "$i".bam -o "$i"-PROFILE
    done
    
    
5. merge

    #!/bin/bash
     anvi-merge *PROFILE/PROFILE.db -c anvi-contigs.db -o 2MPRE-MERGED
     
6. After merging, the contigs were binned manually using Anvio. Files required to load an interactive anvio display that includes the placement of contigs into bins can be found here https://figshare.com/account/home#/projects/77799

### Each of the bins collected from the 27 independetly assembled samples were moved to a single directory and characterized for taxonomy, dereplicated, and annotated. The steps for each of these processes are described below.


