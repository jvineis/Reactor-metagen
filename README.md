# Reactor-metagen
## Documentation of the metagenomic analysis including MAG collections, MAG annotation, phylogenetics, phylogenomics, and pangenomics. The quality filtered reads used as the input for the assembly and all mapping steps are found here (figshare link in progress). The metadata for each of the samples is contained in the supplemental table (S1) in the manuscript (link to mscript).

1.Assemble the reads. You will need all of the metagenomic datasets and the file "quality-filtered-data-list.txt" in a single directory to run the following command. If you have a SLURM management system, please contact me for alternate methods to run this step through an array submission.
     
    #!/bin/bash
    module load megahit/1.0.6
    for i in cat quality-filtered-data-list.txt
    do
       gunzip "$i".gz
       megahit -r "$i" -o "$i"-MEGAHIT
       gzip "$i"
    done
    
2.Map the reads. Run these commands from within the directory containing the megahit assembly. The commands are also looking for the quality filtered reads in a directory above which will also need to contain the txt file "abm-samples.txt"
    
    #!/bin/bash
    for i in `cat ../abm-samples.txt`
    do
        bowtie2 -x anvi-contigs -f -U ../"$i" -S "$i".sam
        samtools view -bS -F 4 "$i".sam > "$i"-RAW.bam
        rm -rf "$i".sam
        anvi-init-bam "$i"-RAW.bam -o "$i".bam
        rm -rf "$i"-RAW.bam
    done

3.Build the Anvio Contigs and the Bowtie2 dbs.

    #!/bin/bash
    anvi-script-reformat-fasta final.contigs.fa -o anvi-contigs.fa -l 2000 --simplify-names
    bowtie2-build anvi-contigs.fa anvi-contigs
    anvi-gen-contigs-database -f anvi-contigs.fa -o anvi-contigs.db
    anvi-run-hmms -c anvi-contigs.db
     
4.profile with Anvio. You will need to run this command separately for each of the assemblies. 

    #!/bin/bash
    for i in `cat ../jhv-samples.txt`
    do
        clusterize anvi-profile -c anvi-contigs.db -T 5 -i "$i".bam -o "$i"-PROFILE
    done  
    
5.Merge the Anvio profiles for each sample independently. This command will need to be run separately for each sample.  

    #!/bin/bash
     anvi-merge *PROFILE/PROFILE.db -c anvi-contigs.db -o 2MPRE-MERGED
     
6.After merging, the contigs were binned manually using Anvio. Files required to load an interactive anvio display that includes the placement of contigs into bins can be found here https://figshare.com/account/home#/projects/77799

### Each of the bins collected from the 27 independetly assembled samples were moved to a single directory and characterized for taxonomy, dereplicated, and annotated. The steps for each of these processes are described below. I have included the SLURM submission details.

1.Taxonomy assignment via GTDBtk: version 0.3.2 Copyright 2017 Pierre Chaumeil. The bactchfile-genomes.txt file is located in this repository. 

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=10
    #SBATCH --time=4:00:00
    #SBATCH --mem=300Gb
    #SBATCH --partition=short

    gtdbtk classify_wf --batchfile batchfile-genomes.txt -x fa --out_dir GTDBtk-OUTPUT
    
### Pangenomics for Chlorobium and Sedimenticola.

1.First you need to search for the Chlorobium genomes that are contained at IMG https://img.jgi.doe.gov/cgi-bin/mer/main.cgi and NCBI (e.g. https://www.ncbi.nlm.nih.gov/assembly/?term=Sedimenticola)
and download them to your server. Depending on where/how you download the references, you will need to decompress the files etc. The main objective of the first step is to get all the genomes you want to include in your pangenomic analysis ino a single directory

For example, the result for the NCBI search for Sedimenticola resulted in a tar archive. Here is how I decompress and move the fasta files to a single directory where I will stage all of the fastas that I want to include in the pangenomic analysis.

    mkdir /work/jennifer.bowen/SEDIMENTICOLA-PAN/GENOME-FASTAS
    rsync -HalP ~/Downloads/genome_assemblies_genome_fasta.tar vineis.j@login.discovery.neu.edu:/work/jennifer.bowen/SEDIMENTICOLA-PAN/GENOME-FASTAS/
    tar -xvf genome_assemblies_genome_fasta.tar
    cd ncbi-genomes-2021-01-12
    mv *.fna.gz ../
    rm -rf ncbi-genomes-2021-01-12
    
2.So, step 1 gave you the power to get the genomes yourself, but the genomes used for the pangeomic analysis can be obtained from here (link to figshare) so that you can conduct the same analysis found in the paper without having to collect them from all over the place like I did. The approach found here is the same for both Chlorobium and Sedimenticola pangenomic analysis. This is everything you need to do in order to characterize the genomes/mags in your directory.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=04:00:00
    #SBATCH --mem=200Gb
    #SBATCH --partition=short
    #SBATCH --array=1-5

    #Anvio- 7.0

    ASSEMBLY=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_jgi-genomes.txt)
    anvi-script-reformat-fasta ${ASSEMBLY}.fna --simplify-names -r jgi_${ASSEMBLY}-simplify-report.txt -o jgi_${ASSEMBLY}.fa
    anvi-gen-contigs-database -f jgi_${ASSEMBLY}.fa -o jgi_${ASSEMBLY}.db
    anvi-run-hmms -c jgi_${ASSEMBLY}.db -T 10
    anvi-run-scg-taxonomy -c jgi_${ASSEMBLY}.db -T 30
    anvi-run-hmms -c jgi_${ASSEMBLY}.db -T 10 -H ~/scripts/databas/HMM_co2fix/
    anvi-run-hmms -c jgi_${ASSEMBLY}.db -T 10 -H ~/scripts/databas/all_fungene_anvio/
    anvi-run-kegg-kofams -c jgi_${ASSEMBLY}.db -T 40
    anvi-run-ncbi-cogs -c jgi_${ASSEMBLY}.db -T 40
    anvi-run-pfams -c jgi_${ASSEMBLY}.db -T 40
    
3.Now you are pretty much ready to run the pangenomic analysis. Just need to build a file called external-genomes.txt like this
run from the same GENOME-FASTAS directory that you have been working from.

    ls *.db | cut -f 1 -d "." > 1
    ls *.db > 2
    paste 1 2 > external-genomes.txt
    
add a header to this file and you are good to go for the next step.  The best looking file ever looks like this with a bunch more genomes and corresponding contigs_db_path below.  

    name    contigs_db_path
    GCA_000428045	GCA_000428045.db
    GCA_001007875	GCA_001007875.db
    GCA_002084045	GCA_002084045.db
    GCA_002084055	GCA_002084055.db

Here is how you run the pangenomic analysis.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=05:00:00
    #SBATCH --mem=200Gb
    #SBATCH --partition=short

    anvi-gen-genome-storage -e external-genomes.txt -o SEDIMENTICOLA-GENOMES.db
    anvi-pan-genome -g SEDIMENTICOLA-GENOMES.db --project-name "Sedimenticola_PAN" --output-dir SEDIMENTICOLA-PAN --num-threads 10 --minbit 0.5 --mcl-inflation 10 --min-occurrence 2 --use-ncbi-blast





