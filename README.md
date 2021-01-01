# Reactor-metagen
## Documentation of the metagenomic analysis including MAG collections, MAG annotation, phylogenetics, phylogenomics, and pangenomics.

## Lets begin with MAG recovery. Here are the general steps used. The location of directories may be invalid for some of the work here due to the removal of temporary files etc.  Its best to start from the beginning if you are not sure.  

1. Assemble the reads.  
     
    #!/bin/bash
    module load megahit/1.0.6
    for i in `cat abm-sulfate-enriched-samples.txt`
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
