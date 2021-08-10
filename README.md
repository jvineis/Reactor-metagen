# Reactor-metagen
## Documentation of the metagenomic analysis including MAG collections, MAG annotation, phylogenetics, phylogenomics, and pangenomics. The quality filtered reads used as the input for the assembly and all mapping steps are found at MGRAST mgp84173 (https://www.mg-rast.org/mgmain.html?mgpage=search&search=mgp84173). The metadata for each of the samples is contained in this git (Table-S1.xlsx).

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
    
 2.Dereplication
     
    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=10
    #SBATCH --time=24:00:00
    #SBATCH --mem=300Gb
    #SBATCH --partition=short

    anvi-compute-genome-similarity -e external-genomes.txt -o x_ANI -T 40
    anvi-dereplicate-genomes --ani-dir x_ANI/ -o x_ANI_dereplication --program pyANI --method ANIb --use-full-percent-identity --min-full-percent-identity 0.90 --similarity-threshold 0.95
    
3.Run DRAM

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=10
    #SBATCH --time=24:00:00
    #SBATCH --mem=250Gb
    #SBATCH --partition=short

    DRAM.py annotate -i '*fa' -o x_DRAM-annotate --gtdb_taxonomy gtdbtk.bac120.summary.tsv

4.Run fungene HMMs

### Relative abundance calculations for each MAG

1.First you will need a list of the dereplicated MAGs so that you only map to one representative MAG. Using more than one of the redundant genomes that represents the same organism in a sample will cause inaccurate coverage and depth estimates. So run this to get the list of unique MAGs.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=10
    #SBATCH --time=24:00:00
    #SBATCH --mem=300Gb
    #SBATCH --partition=short


    anvi-dereplicate-genomes --ani-dir x_ANI/ -o x_ANI_dereplication --program pyANI --method ANIb --use-full-percent-identity --min-full-percent-identity 0.90 --similarity-threshold 0.95

2.Pull out the list of representative MAGs from the cluster report

    cut -f 3 x_ANI_dereplication/CLUSTER_REPORT.txt > x_DEREPLICATED-MAG-IDs.txt
      
3.Fix the deflines of the fasta to reflect the name of the MAG within each of the deflines using the x_append-name-to-fasta-deflines.py script

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=2
    #SBATCH --time=1:00:00
    #SBATCH --mem=200Gb
    #SBATCH --partition=express
    #SBATCH --array=1-113

    binfa=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_samples.txt)
    python x_append-name-to-fasta-deflines.py ${binfa}

3.copy all of the mags in the unique list to the mapping directory.

    for i in `cat x_DEREPLICATED-MAG-IDs.txt`; do cp $i'-bdeflines.fa' MAPPING; done
    cd MAPPING
    cat *.fa > x_ALL-DEREP-MAGS.fa
    rm *bdeflines.fa
    
4.Run the mapping script using bbmap

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=10
    #SBATCH --time=22:00:00
    #SBATCH --mem=300Gb
    #SBATCH --partition=short
    #SBATCH --array=1-33

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_sample-names.txt)

    bbmap.sh ref=x_ALL-DEREP-MAGS.fa nodisk=true in=/work/jennifer.bowen/JOE/FTR/QUALITY-FILTERED-READS/${SAMPLE} covstats=MAPPING/${SAMPLE}-covstats.txt scafstats=MAPPING/${SAMPLE}-scafstats.txt out=MAPPING/${SAMPLE}.bam bamscript=to_bam.sh

5.collect a list of the scaffold ids.. required for tabulating coverage of each scaffold.

    grep ">" x_ALL-DEREP-MAGS.fa | sed 's/>//g' > x_ALL-DEREP-MAGS-scaf-ids.txt

6.Tabulate the total number of reads mapped to the MAG using estimate-MAG-coverage-from-bbmap-covstats-v2-FTR.py

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=2
    #SBATCH --time=1:00:00
    #SBATCH --mem=20Gb
    #SBATCH --partition=express
    #SBATCH --array=1-33

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_sample-names.txt)
    python ~/scripts/estimate-MAG-coverage-from-bbmap-covstats-v2-FTR.py --mappingfile ${SAMPLE}-scafstats.txt --out ${SAMPLE}-bases-recruited-per-MAG.txt --ids x_ALL-DEREP-MAGS-scaf-ids.txt
 
7.Paste together the individual results of the mapping, open in excel and go to work calculating the normalized relative abundance using the number of bases in the MAG (length) and the total number of bases for each sample.

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

4.So you are now ready to compute ANI among the MAGs.

    #!/bin/bash
    #    
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=05:00:00
    #SBATCH --mem=200Gb
    #SBATCH --partition=short

    anvi-compute-genome-similarity -e external-genomes.txt -o ANI_refined -p SEDIMENTICOLA-PAN/Sedimenticola_PAN-PAN.db -T 30

5.Add some data to the pangenomic db

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=50Gb
    #SBATCH --partition=express

    anvi-import-misc-data Chlorobium-samples-metadata-for-anvio.txt -p CHLOROBIUM-PAN/Chlorobium_PAN-PAN.db -t layers

6.To display the pangenome, simply run this to log into your server

    anvi-display-pan -g SEDIMENTICOLA-refined-GENOMES.db -p SEDIMENTICOLA-PAN-refined/Sedimenticola_PAN-PAN.db --server-only -P 8087


Phylogenomics of FTR MAGs and Genomes of Earth's Microbiomes (GEM)s



