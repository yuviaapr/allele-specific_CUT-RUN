# Allele-specific CUT&RUN workflow

Nextflow workflow to process raw TF or histone CUT&RUN data with control libraries for allele-specific analysis. It takes compressed fastq files as input and generates quality metrics, normalized signal tracks, peak files based on total signal ("Gall") using macs2 and/or SEACR, and bed files with nucleosome and nucleosome-free fragments. Also, it saves log files that can be parsed to calculate the number of reads filtered out in each step.

## Author

Yuvia Alhelí PÉREZ-RICO

Affiliation: European Molecular Biology Laboratory | EMBL

## Dependencies

To use this workflow, you will need to install the following programs, the indicated version is the one that I have used to analyze my data:

- nextflow (20.04.1)
- FastQC (v0.11.9)
- Trim Galore (0.6.3)
- Bowtie2 (2.3.4.1)
- samtools (1.9)
- SNPsplit (0.3.4)
- Picard Tools (2.20.8)
- deepTools (3.4.3)
- bedtools (v2.29.2)
- macs2 (2.2.7.1)
- SEACR (1.3)
- R (4.0.2)
- ATACseqQC bioconductor package (1.12.0)
- bioconductor packages for peak annotation: ggupset (0.3.0), org.Mm.eg.db (3.11.4), ChIPseeker (1.24.0) and TxDb.Mmusculus.UCSC.mm10.knownGene (3.10.0).
- bedGraphToBigWig (v 4)
- featureCounts (v2.0.1)
- R scripts included as part of my repository for [ChIP-seq](https://github.com/yuviaapr/allele-specific_ChIP-seq/tree/main/bin) data processing: plot_fragment_size_distribution.R and plot_peak_distribution_features.R

Note: I suggest to install all programs in a conda environment, except the packages used for peak annotation, which I installed in an independent environment and call it within the nextflow script because I experienced dependency issues (line 952 of main script).

## How to use the workflow

Thank you for your interest in using the workflow!

After downloading the main script, the configuration file and additional scripts, you will need to prepare the following files and change the paths in the configuration file accordingly:

- N-masked bowtie2 index using SNPs of 2 strains.
- Compressed file with the SNPs of the 2 strains for SNPsplit.
- Bed file containing blacklist genomic regions.
- Tab-separated file with chromosome sizes (only 2 columns: chromosome name and size).
- Comma-separated file indicating the sample origin, name, paths to the fastq files and the library type (TF, histone or control). An example is available in the 'docs' folder. The format of the first two columns of this file is important to properly pair target and control libraries, and it is required to include 1 underscore in the 'clone' column.

If you are not using SLURM, then do additional modifications to the configuration file considering the workload manager that you use. Also, indicate if you are analyzing TF data, histone data or both. Finally, write a simple bash script that will be submitted to the cluster to activate the conda environment and start the main nextflow job:

`source /home/user/miniconda2/bin/activate /home/user/conda-envs/CaR_bioc`

`nextflow run allelic_CUTandRUN.nf`

Please, keep in mind that this workflow was written for the processing of mouse paired-end data, therefore, you will need to do some changes to the main and additional scripts before using it with data of other species, for example, changing the strain names and databases for peak annotation. Also, I skip some chromosomes for the normalization of signal tracks, as I know that they tend to get duplicated in the cell lines that are commonly handled in our lab or to ease the visualization of differentiated and undifferentiated cells. Modify the 'nucleosome_tracks' and 'NFR_tracks'  processes in case that you do not want to skip chromosomes.

As I am interested in obtaining allele-specific signal, the CUT&RUN libraries that I generally process are deeply sequenced (> 30M) compared to standard CUT&RUN libraries. I noticed that when using all fragments to call peaks with SEACR, the resulting peaks tend to be very large and expand regions with low signal. For this reason, I downsample histone data to 10 M fragments using the 'numFrag' parameter that can be modified in the configuration file.

## Acknowledgements

This repository is part of a project that has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 882771 – TAXi-ch.

## License

Licenced under the GNU general public license.

