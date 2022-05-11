#!/usr/bin/env nextflow

/*
 Pipeline to process allele-specific CUT&RUN data and generate processed bam and bed files, peaks, read counts and normalized signal tracks
 Author:
	- Yuvia A. PEREZ RICO <yuvia.perez-rico@embl.de>
*/

log.info "     allelic CUTandRun processing - version 0.1.0      "
log.info "#######################################################"
log.info "Sample description file	= ${params.sampleInfo}"
log.info "Reads per split file		= ${params.chunkSize}"
log.info "Mouse strain 1		= ${params.G1}"
log.info "Mouse strain 2		= ${params.G2}"
log.info "Output directory		= ${params.outDir}"
log.info "Path to genome index		= ${params.genomeDirPath}"
log.info "SNP reference file		= ${params.snpFile}"
log.info "Blacklisted regions		= ${params.blackList}"
log.info "Script to plot fragment sizes	= ${params.rsFrag}"
log.info "Script to annotate peaks	= ${params.rsPeaks}"
log.info "TF data included		= ${params.dataTF}"
log.info "Histone data included		= ${params.dataHistone}"
log.info "SEACR script			= ${params.shSEACR}"
log.info "Chromosome sizes		= ${params.chrmSizes}"
log.info "Fragments down-sampling	= ${params.numFrag}"
log.info "Number of threads		= ${params.numCPUs}"
log.info "Number of threads deepTools	= ${params.numCPUs_Dtools}"
log.info "\n"

/*
 Validate input parameters
*/

if( !(params.chunkSize instanceof Number) ){
	exit 1, "Invalid chunk size = ${params.chunkSize}"
}

if( !(params.numFrag instanceof Number) ){
	exit 1, "Invalid number of fragments for down-sampling = ${params.numFrag}"
}

if( !(params.numCPUs instanceof Number) ){
	exit 1, "Invalid number of CPUs = ${params.numCPUs}"
}

if( !(params.numCPUs_Dtools instanceof Number) ){
	exit 1, "Invalid number of CPUs for deepTools = ${params.numCPUs_Dtools}"
}

if( !(params.G1 in ['C57BL-6J', 'CAST-EiJ', 'PWK-PhJ', '129S1-SvImJ', 'FVB-NJ', '129P2-OlaHsd']) ){
	exit 1, "Invalid strain 1 name = ${params.G1}"
}

if( !(params.G2 in ['C57BL-6J', 'CAST-EiJ', 'PWK-PhJ', '129S1-SvImJ', 'FVB-NJ', '129P2-OlaHsd']) ){
	exit 1, "Invalid strain 2 name = ${params.G2}"
}

if( !(params.dataTF in ['yes', 'no']) ){
	exit 1, "Invalid option for TF data = ${params.dataTF}"
}

if( !(params.dataHistone in ['yes', 'no']) ){
	exit 1, "Invalid option for histone data = ${params.dataHistone}"
}

/*
 Validate input files
*/

sdFile = file(params.sampleInfo)
if( !sdFile.exists() ){
	exit 1, "The specified sample description file does not exist = ${params.sampleInfo}"
}
log.info "Checking sample description file = $sdFile"

varFile = file(params.snpFile)
if( !varFile.exists() ){
	exit 1, "The specified SNP annotation file does not exist = ${params.snpFile}"
}
log.info "Checking SNP annotations file = $varFile"

chrFile = file(params.chrmSizes)
if( !chrFile.exists() ){
	exit 1, "The specified file with the chromosome sizes does not exist = ${params.chrmSizes}"
}
log.info "Checking chromosome information = $chrFile"

blReg = file(params.blackList)
if( !blReg.exists() ){
	exit 1, "The specified blacklist file does not exist = ${params.blackList}"
}
log.info "Checking blackList regions = $blReg"

plotFGM = file(params.rsFrag)
if( !plotFGM.exists() ){
	exit 1, "The specified R script does not exist = ${params.rsFrag}"
}
log.info "Checking R script to generate fragment size distribution = $plotFGM"

plotPeakAnno = file(params.rsPeaks)
if( !plotPeakAnno.exists() ){
	exit 1, "The specified R script does not exist = ${params.rsPeaks}"
}
log.info "Checking R script to annotate peaks = $plotPeakAnno"

SEACR = file(params.shSEACR)
if( !SEACR.exists() ){
	exit 1, "The specified script to run SEACR does not exist = ${params.shSEACR}"
}
log.info "Checking SEACR script to call histone mark regions = $SEACR"

resDir = file(params.outDir)
if( !resDir.exists() && !resDir.mkdirs() ){
	exit 1, "The specified directory to save results cannot be created = ${params.outDir}\n Check file system access permission"
}
log.info "Checking results directory = $resDir"

/*
 Program versions
*/

process get_program_versions{
	publishDir "${resDir}/software", mode: 'move'

	output:
	file('programs_version.txt') into programs_version

	"""
	echo nextflow ${nextflow.version} > tmp_version.txt
	fastqc -v >> tmp_version.txt
	echo trim_galore \$(trim_galore -v | awk '/version/{print\$2}') >> tmp_version.txt
	echo bowtie2 \$(bowtie2 --version | head -1 | cut -d " " -f 3) >> tmp_version.txt
	samtools --version | grep samtools >> tmp_version.txt
	echo SNPsplit \$(SNPsplit --version | awk '/Version/{print\$2}') >> tmp_version.txt
	picard SortSam --version 2>&1 | awk '{sub(/-SNAPSHOT/,"");print"picard "\$1}' >> tmp_version.txt
	bamCoverage --version >> tmp_version.txt
	bedtools --version >> tmp_version.txt
	macs2 --version >> tmp_version.txt
	R --version | head -1 | cut -f 1-3 -d " " >> tmp_version.txt
	conda list -f bioconductor-atacseqqc | grep bioconda | awk '{print\$1" "\$2}' >> tmp_version.txt
	echo ${SEACR} | awk '{sub(/_/," ");sub(/\\.sh/,"");print"SEACR "\$2" "\$3}' >> tmp_version.txt
	bedGraphToBigWig 2>&1 | grep Convert | cut -f 1 -d "-" >> tmp_version.txt
	featureCounts -v 2>&1 | grep feature >> tmp_version.txt
	sort tmp_version.txt > programs_version.txt
	"""

}

/*
 Create channels with the fastq files for processing and quality control
*/

// Reads1

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.clone + "_" + row.sample, file(row.reads1)) }
	.set { samples }

samples.into { samples_r1; samples_count }

// Reads2

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.clone + "_" + row.sample, file(row.reads2)) }
	.set { samples_r2 }

// Data type

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.clone, row.sample, row.type) }
	.into { dataType1; dataType2; dataType3 }

// Both reads

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.clone + "_" + row.sample, file(row.reads1), file(row.reads2)) }
	.set { samples_quality }

/*
 Step 0. Pre-alignment quality check
*/

process fastq_quality_pre{
	publishDir "${resDir}/qc/1_fastqc_pre-alignment", mode: 'copy'
	label 'fastqQual'

	input:
	set val(name), file(reads1), file(reads2) from samples_quality

	output:
	file("${name}_fastqc") into fastqc_pre_results

	"""
	mkdir ${name}_fastqc
	fastqc -o ${name}_fastqc -q -t ${params.numCPUs} ${reads1} ${reads2}
	"""

}

/*
 Step 1. Count raw number of read pairs
*/

process raw_pairs{
	input:
	set val(name), file(reads) from samples_count

	output:
	set val(name), stdout into total_pairs

	"""
	echo \$(zcat ${reads} | wc -l)/4 | bc | tr -dc '0-9'
	"""

}

/*
 Step 2. Split fastq files to ease setting the configuration file
*/

process split_Reads1{
	label 'fastq_splitting'

	input:
	set val(name), file(reads) from samples_r1

	output:
	set val(name), file('*.gz') into samples_r1_split

	"""
	zcat ${reads} | split --numeric-suffixes=1 -a 4 -l \$((${params.chunkSize} * 4)) --filter='gzip > \$FILE.gz' - "${name}_reads1_"
	"""

}

process split_Reads2{
	label 'fastq_splitting'

	input:
	set val(name), file(reads) from samples_r2

	output:
	set val(name), file('*.gz') into samples_r2_split

	"""
	zcat ${reads} | split --numeric-suffixes=1 -a 4 -l \$((${params.chunkSize} * 4)) --filter='gzip > \$FILE.gz' - "${name}_reads2_"
	"""

}

// Organize paired-end reads

samples_r1_split
	.combine(samples_r2_split, by: 0)
	.transpose()
	.set { samples }

/*
 Step 3. Trim adapters and low quality bases
*/

process trim_reads{
	publishDir "${resDir}/qc/2_trimming", mode: 'copy', pattern: "${name}/*report.txt"

	input:
	set val(name), file(reads1), file(reads2) from samples

	output:
	set val(name), file("${name}/*{1,2}.fq.gz") into trimmed_samples
	file("${name}/*report.txt") into trim_stats

	"""
	trim_galore --stringency 1 --length 25 --quality 10 --paired ${reads1} ${reads2} -o ${name}
	"""

}

/*
 Step 4. Mapping
*/

process read_mapping{
	publishDir "$resDir/qc/3_mapping", mode: 'copy', pattern: '*.log'

	input:
	set val(name), file(trimmed_reads) from trimmed_samples

	output:
	set val(name), file('*.bam') into split_mapping
	file('*.log') into mapping_stats

	"""
	file=\$(echo ${trimmed_reads[0]} | sed s/.gz_trimmed.fq.gz//)
	bowtie2 --very-sensitive -X 1000 --no-mixed --no-discordant -p ${params.numCPUs} -x ${params.genomeDirPath} \
		-1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} 2> \$file.log | samtools view -b > \$file.bam
	"""

}

/*
 Step 5. Select properly paired reads with high mapping quality
*/

process pairs_filter{
	publishDir "${resDir}/qc/4_pairs_filter", mode: 'copy', pattern: '*.stats'
	label 'filter_bams'

	input:
	set val(name), file(mapped_reads) from split_mapping

	output:
	set val(name), file('*.bam') into split_mapping_pairs
	file('*.stats') into pairs_stats

	"""
	samtools flagstat -@ ${params.numCPUs} ${mapped_reads} > ${mapped_reads.baseName}_pre-rmQC.stats
	samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 -q 20 -@ ${params.numCPUs} ${mapped_reads} > ${mapped_reads.baseName}_rmQC.bam
	"""

}

/*
 Step 6. Remove mitochondrial reads
*/

process chrM_filter{
	publishDir "${resDir}/qc/4_pairs_filter", mode: 'copy', pattern: '*QC.stats'
	publishDir "${resDir}/qc/5_chrM_filter", mode: 'copy', pattern: '*ChrM.stats'
	label 'filter_bams'

	input:
	set val(name), file(mapped_pairs) from split_mapping_pairs

	output:
	set val(name), file('*.bam') into split_mapping_noMit, split_mapping_noMit_pbc
	file('*QC.stats') into postQC_stats
	file('*ChrM.stats') into chrM_stats

	"""
	samtools flagstat -@ ${params.numCPUs} ${mapped_pairs} > ${mapped_pairs.baseName}_post-rmQC.stats
	samtools view -h -@ ${params.numCPUs} ${mapped_pairs} | grep -v chrM | samtools view -b -@ ${params.numCPUs} - > ${mapped_pairs.baseName}_rmChrM.bam
	samtools flagstat -@ ${params.numCPUs} ${mapped_pairs.baseName}_rmChrM.bam > ${mapped_pairs.baseName}_rmChrM.stats
	"""

}

/*
 Step 7. Calculate PCR bottlenecking coefficients 1 and 2
*/

// Organise files

split_mapping_noMit_pbc.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { bam_pbc }

process PCR_coeffs{
	publishDir "${resDir}/qc/6_PCR_bottlenecking_coefficients", mode: 'copy', pattern: '*PBC2.txt'

	input:
	set val(name), file(filtered_pairs) from bam_pbc

	output:
	file('*PBC2.txt') into PBCs

	"""
	samtools merge -f ${name}.bam -@ ${params.numCPUs} ${filtered_pairs}
	samtools sort -n -O BAM -o ${name}_sortName.bam -@ ${params.numCPUs} ${name}.bam
	bedtools bamtobed -bedpe -i ${name}_sortName.bam | awk '{print\$1"\t"\$2"\t"\$4"\t"\$6"\t"\$9"\t"\$10}' | sort | uniq -c | \
		awk 'BEGIN{M1=0;M2=0;MD=0;MT=0;print"M1\tM2\tMdistinct\tMtotal\tPBC1\tPBC2"}; \$1==1 {M1=M1+1}; \
		\$1==2 {M2=M2+1}; {MD=MD+1;MT=MT+\$1}; END{print M1"\t"M2"\t"MD"\t"MT"\t"M1/MD"\t"M1/M2}' > ${name}_PBC1_PBC2.txt
	"""

}

/*
 Step 8. Genotype assignation
*/

process SNPsplit{
	publishDir "${resDir}/qc/7_SNPsplit", mode: 'copy', pattern: "${name}/*.txt"

	input:
	set val(name), file(filtered_pairs) from split_mapping_noMit

	output:
	set val(name), file("${name}/*genome1.bam") into split_genome1
	set val(name), file("${name}/*genome2.bam") into split_genome2
	set val(name), file("${name}/*unassigned.bam") into split_unassigned
	file("${name}/*report.txt") into SNPs_report_stats
	file("${name}/*sort.txt") into SNPs_sort_stats

	"""
	SNPsplit --paired --no_sort --snp_file ${varFile} -o ${name} ${filtered_pairs}
	"""

}

/*
 Step 9. Merge and sort split BAM files per sample and remove duplicates from merged BAM file
*/

// Organise files

split_genome1.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { genome1 }

split_genome2.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { genome2 }

split_unassigned.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { unassigned }

// Process G1 files

process filter_G1_bams{
	publishDir "${resDir}/processed_bams", mode: 'copy', pattern: '*rmdup.bam'
	publishDir "${resDir}/qc/8_duplicates", mode: 'copy', pattern: '*rmdup.stats'
	label 'process_bams'

	input:
	set val(name), file(genome1_bam_files) from genome1

	output:
	set val(name), file('*rmdup.bam') into G1_bam_merge, G1_bam, G1_count, G1_coverage
	file('*rmdup.stats') into duplicates_stats_G1

	"""
	samtools merge -f ${name}_${params.G1}.bam ${genome1_bam_files}
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_${params.G1}.bam O=${name}_${params.G1}_sorted.bam SORT_ORDER=coordinate
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MarkDuplicates I=${name}_${params.G1}_sorted.bam O=${name}_${params.G1}_sorted_rmdup.bam \
		 M=${name}_${params.G1}_sorted_rmdup.stats REMOVE_DUPLICATES=true ASSUME_SORTED=true
	"""

}

// Process G2 files

process filter_G2_bams{
	publishDir "${resDir}/processed_bams", mode: 'copy', pattern: '*rmdup.bam'
	publishDir "${resDir}/qc/8_duplicates", mode: 'copy', pattern: '*rmdup.stats'
	label 'process_bams'

	input:
	set val(name), file(genome2_bam_files) from genome2

	output:
	set val(name), file('*rmdup.bam') into G2_bam_merge, G2_bam, G2_count, G2_coverage
	file('*rmdup.stats') into duplicates_stats_G2

	"""
	samtools merge -f ${name}_${params.G2}.bam ${genome2_bam_files}
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_${params.G2}.bam O=${name}_${params.G2}_sorted.bam SORT_ORDER=coordinate
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MarkDuplicates I=${name}_${params.G2}_sorted.bam O=${name}_${params.G2}_sorted_rmdup.bam \
		 M=${name}_${params.G2}_sorted_rmdup.stats REMOVE_DUPLICATES=true ASSUME_SORTED=true
	"""

}

// Process unassigned files

process filter_unassigned_bams{
	publishDir "${resDir}/qc/8_duplicates", mode: 'copy', pattern: '*rmdup.stats'

	input:
	set val(name), file(unassigned_bam_files) from unassigned

	output:
	set val(name), file('*rmdup.bam') into unassigned_bam_merge
	file('*rmdup.stats') into duplicates_stats_unassigned

	"""
	samtools merge -f ${name}_unassigned.bam ${unassigned_bam_files}
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_unassigned.bam O=${name}_unassigned_sorted.bam SORT_ORDER=coordinate
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MarkDuplicates I=${name}_unassigned_sorted.bam O=${name}_unassigned_sorted_rmdup.bam \
		 M=${name}_unassigned_sorted_rmdup.stats REMOVE_DUPLICATES=true ASSUME_SORTED=true
	"""

}

// Merge all read pairs

process merge_mapped{
	publishDir "${resDir}/processed_bams", mode: 'copy', pattern: '*.bam'
	publishDir "${resDir}/qc/9_processed_bams", mode: 'copy', pattern: '*.stats'

	input:
	set val(name), file(genome1_bam), file(genome2_bam), file(unassigned_bam) from G1_bam_merge.join(G2_bam_merge).join(unassigned_bam_merge)

	output:
	set val(name), file('*Gall_sorted.bam') into bam_fraction, all_bam, bam_quality, bam_peaks
	file('*.stats') into Gall_stats

	"""
	samtools merge -f ${name}_Gall.bam ${genome1_bam} ${genome2_bam} ${unassigned_bam}
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_Gall.bam O=${name}_Gall_sorted.bam \
		SORT_ORDER=coordinate
	samtools flagstat ${name}_Gall_sorted.bam > ${name}_processed_bam.stats
	"""

}

// Set combined channel

G1_bam.mix(G2_bam).mix(all_bam).into { bam_coverage; bam_read_counts }

// Replicate channel for quality checks

bam_quality.into { bam_post; bam_fgm; bam_count }

/*
 Step 10. Calculate the percentage of reads assigned to each allele
*/

process allelic_reads{
	publishDir "${resDir}/qc/10_allele-specific_reads", mode: 'copy'
	label 'stats'

	input:
	set val(name), file(genome1_bam), file(genome2_bam), file(Gall_bam) from G1_count.join(G2_count).join(bam_count)

	output:
	file('*reads.txt') into allelic_reads

	"""
	samtools view -c ${genome1_bam} > ${name}_reads_per_allele.txt
	samtools view -c ${genome2_bam} >> ${name}_reads_per_allele.txt
	total=\$(samtools view -c ${Gall_bam})
	awk -v all=\$total 'BEGIN{print"allele\tassigned reads\tpercentage";c=1};{print"G"c"\t"\$1"\t"(\$1*100)/all;c=c+1}; \
		END{print"Total reads\t"all}' ${name}_reads_per_allele.txt > ${name}_allelic_reads.txt
	"""

}

/*
 Step 11. Post-alignment quality check
*/

process fastq_quality_post{
	publishDir "${resDir}/qc/11_fastqc_post-alignment", mode: 'copy'
	label 'fastqQual'

	input:
	set val(name), file(Gall_bam) from bam_post

	output:
	file("${name}_fastqc") into fastqc_post_results

	"""
	mkdir ${name}_fastqc
	fastqc -o ${name}_fastqc -q -t ${params.numCPUs} ${Gall_bam}
	"""

}

/*
 Step 12. Plot fragment size distribution
*/

process plot_fragments{
	publishDir "${resDir}/qc/12_fragments_distribution", mode: 'copy'
	label 'scR'

	input:
	set val(name), file(Gall_bam) from bam_fgm

	output:
	file('*.pdf') into plots_fgm

	"""
	samtools index ${Gall_bam}
	Rscript ${plotFGM} ${Gall_bam} ${name}
	"""

}

/*
 Step 13. Calculate the non-redundant fraction
*/

process non_redundant{
	publishDir "${resDir}/qc/13_non-redundant_fraction", mode: 'copy'
	label 'stats'

	input:
	set val(name), file(bam), val(total) from bam_fraction.combine(total_pairs, by: 0)

	output:
	file('*fraction.txt') into NRF

	"""
	samtools flagstat ${bam} | grep read1 | awk -v all=${total} '{print\$1/all}' > ${name}_nonRedundant_fraction.txt
	"""

}

/*
 Step 14. Call peaks using macs2 and SEACR
*/

// Generate channel with matched TF and control data

bam_peaks
	.map{ row -> tuple(row[0].split("_")[0] + "_" + row[0].split("_")[1], row[0].split("_")[2], row[1]) }
	.join(dataType1, by: [0,1])
	.branch {
		control: it[3] == 'control'
		TF: it[3] == 'TF'
		histone: it[3] == 'histone'
	}
	.set { bam }

bam.control.into { bam_control1; bam_control2 }

bam.TF.into { bam_TF1; bam_TF2 }

// Generate channel with TF and control data

bam_TF1
	.mix(bam_control1)
	.into { bam_TF_control; bam_NFR_signal }

// Generate channel with histone, TF and control data

bam.histone
	.mix(bam_control2)
	.mix(bam_TF2)
	.set { bam_process }

// Process bam files to extract nucleosomal fragments

process bam_to_bed{
	publishDir "${resDir}/processed_beds", mode: 'copy'
	label 'bam_fragments'

	input:
	set val(clone), val(target), file(bam), val(type) from bam_process

	output:
	set val(clone), val(target), val(type), file("*.bed") into bed_peaks
	set val(clone), val(target), stdout into size_factors_nucleosome

	"""
	samtools sort -@ ${params.numCPUs} -n ${bam} -o ${bam.baseName}Name.bam
	bedtools bamtobed -bedpe -i ${bam.baseName}Name.bam | awk '\$1==\$4 && \$6-\$2 > 120 && \$6-\$2 <= 500 {print\$1"\t"\$2"\t"\$6}' | \
		sort -k1,1 -k2,2n -k3,3n > ${bam.baseName}_nucleosome.bed
	echo 100000000/\$(wc -l ${bam.baseName}_nucleosome.bed | cut -d " " -f 1) | bc -l
	"""

}

// Combine histone and TF data with their respective controls

bed_peaks
	.branch {
		control: it[2] == 'control'
		TF: it[2] == 'TF'
		histone: it[2] == 'histone'
	}
	.set { beds }

beds.control.into { bed_control1; bed_control2; bed_control3 }

beds.TF.into { bed_TF1; bed_TF2 }

beds.histone.into { bed_histone1; bed_histone2; bed_histone3 }

// Generate channel with histone and control data for macs2

bed_histone1
	.combine(bed_control1, by: 0)
	.map{ row -> tuple(row[0], row[1], row[2], row[3], row[4], row[6]) }
	.set { bed_histone_macs }

// Generate channel with TF and control data for macs2

bed_TF1
	.combine(bed_control2, by: 0)
	.map{ row -> tuple(row[0], row[1], row[2], row[3], row[4], row[6]) }
	.into { bed_TF_macs; bed_TF_tracks }

// Generate channel with histone and control data for SEACR

bed_histone2
	.mix(bed_control3)
	.set { bed_histone_seacr }

// Generate channel with TF and histone data

bed_TF2
	.mix(bed_histone3)
	.map{ row -> tuple(row[0], row[1], row[3]) }
	.set { bed_frip_macs }

// Call peaks with macs2

process peaks_macs2{
	publishDir "${resDir}/peak_calling", mode: 'copy', pattern: '*_macs2'
	label 'peak_calling'

	input:
	set val(clone), val(target), val(type), file(bed_target), val(control), file(bed_control) from bed_TF_macs.mix(bed_histone_macs)

	output:
	file("${clone}_${target}_${type}_macs2") into macs2_results
	set val(clone), val(target), val(p_caller), file("${clone}_${target}_${type}_macs2/*rmQval.*Peak") into peaks_annot_macs, peaks_macs
	set val(clone), val(target), val(control), val(p_caller), file("${clone}_${target}_${type}_macs2/*rmQval.SAF") into peaks_counting_macs

	script:
	p_caller = "macs2"
	if(type == "TF") {
		"""
		macs2 callpeak -f BEDPE -g mm --keep-dup all --outdir ${clone}_${target}_${type}_macs2 --call-summits \
			-n ${clone}_${target}_${control}_control -t ${bed_target} -c ${bed_control}
		bedtools intersect -v -a ${clone}_${target}_${type}_macs2/*.narrowPeak -b ${blReg} | \
			awk '\$9>=4' > ${clone}_${target}_${type}_macs2/${clone}_${target}_${control}_control_peaks_rmBL_rmQval.narrowPeak
		bedtools intersect -v -a ${clone}_${target}_${type}_macs2/*.bed -b ${blReg} | \
			awk '\$5>=4' > ${clone}_${target}_${type}_macs2/${clone}_${target}_${control}_control_summits_rmBL_rmQval.bed
		cut -f 1-3 ${clone}_${target}_${type}_macs2/${clone}_${target}_${control}_control_peaks_rmBL_rmQval.narrowPeak | uniq | \
			awk -v name=${clone}_${target}_${type} 'BEGIN{count=1;print"GeneID\tChr\tStart\tEnd\tStrand"}; \
			{print name"_"count"\t"\$1"\t"\$2+1"\t"\$3"\t.";count=count+1}' > \
			${clone}_${target}_${type}_macs2/${clone}_${target}_${control}_control_peaks_rmBL_rmQval.SAF
		"""
	} else if(type == "histone") {
		"""
		macs2 callpeak -f BEDPE -g mm --keep-dup all --outdir ${clone}_${target}_${type}_macs2 --broad --broad-cutoff 0.01 \
			-n ${clone}_${target}_${control}_control -t ${bed_target} -c ${bed_control}
		bedtools intersect -v -f 0.50 -a ${clone}_${target}_${type}_macs2/*.broadPeak -b ${blReg} | \
			awk '\$7>=3 && \$9>=3' > ${clone}_${target}_${type}_macs2/${clone}_${target}_${control}_control_peaks_rmBL_rmQval.broadPeak
		cut -f 1-3 ${clone}_${target}_${type}_macs2/${clone}_${target}_${control}_control_peaks_rmBL_rmQval.broadPeak | uniq | \
			awk -v name=${clone}_${target}_${type} 'BEGIN{count=1;print"GeneID\tChr\tStart\tEnd\tStrand"}; \
			{print name"_"count"\t"\$1"\t"\$2+1"\t"\$3"\t.";count=count+1}' > \
			${clone}_${target}_${type}_macs2/${clone}_${target}_${control}_control_peaks_rmBL_rmQval.SAF
		"""
	}

}

// Generate bedgraph files for SEACR

process bed_to_bg{
        publishDir "${resDir}/processed_beds", mode: 'copy', pattern: "*${params.numFrag}.bed"
        publishDir "${resDir}/processed_beds/signal", mode: 'copy', pattern: '*forSEACR.bedgraph'

	input:
	set val(clone), val(target), val(type), file(bed) from bed_histone_seacr

	when:
	params.dataHistone == "yes"

	output:
	file("*${params.numFrag}.bed") into down_sampling
	set val(clone), val(target), val(type), file("*forSEACR.bedgraph") into bg_tracks
	set val(clone), val(target), file("*${params.numFrag}.bed") into bed_histone

	"""
	shuf -n ${params.numFrag} ${bed} | sort -k1,1 -k2,2n > ${bed.baseName}_${params.numFrag}.bed
	bedtools genomecov -bg -i ${bed.baseName}_${params.numFrag}.bed -g ${chrFile} > ${bed.baseName}_${params.numFrag}.bedgraph
	sort -k1,1 -k2,2n ${bed.baseName}_${params.numFrag}.bedgraph > ${bed.baseName}_${params.numFrag}_forSEACR.bedgraph
	"""

}

// Combine histone and TF data with their respective controls

bg_tracks
	.branch {
		control: it[2] == 'control'
		histone: it[2] == 'histone'
	}
	.set { bgs }

bgs.histone
	.combine(bgs.control, by: 0)
	.map{ row -> tuple(row[0], row[1], row[2], row[3], row[4], row[6]) }
	.set { bg_histone_seacr }

// Call peaks with SEACR

process peaks_seacr{
	publishDir "${resDir}/peak_calling/SEACR", mode: 'copy', pattern: '*.bed'
	label 'peak_calling'

	input:
	set val(clone), val(target), val(type), file(bg_target), val(control), file(bg_control) from bg_histone_seacr

	when:
	params.dataHistone == "yes"

	output:
	file("*.bed") into seacr_results
	set val(clone), val(target), val(p_caller), file("*relaxed_rmBL.bed") into peaks_annot_seacr, peaks_histone_seacr
	set val(clone), val(target), val(control), val(p_caller), file("*relaxed_rmBL.SAF") into peaks_counting_seacr

	script:
	p_caller = "SEACR"
	"""
	bash ${SEACR} ${bg_target} ${bg_control} norm relaxed ${clone}_${target}_${type}_${control}_control_SEACR
	bash ${SEACR} ${bg_target} ${bg_control} norm stringent ${clone}_${target}_${type}_${control}_control_SEACR
	bedtools intersect -v -f 0.50 -a ${clone}_${target}_${type}_${control}_control_SEACR.relaxed.bed \
		-b ${blReg} > ${clone}_${target}_${type}_${control}_control_SEACR.relaxed_rmBL.bed
	bedtools intersect -v -f 0.50 -a ${clone}_${target}_${type}_${control}_control_SEACR.stringent.bed \
		-b ${blReg} > ${clone}_${target}_${type}_${control}_control_SEACR.stringent_rmBL.bed
	awk 'BEGIN{print"GeneID\tChr\tStart\tEnd\tStrand"};{print\$6"\t"\$1"\t"\$2+1"\t"\$3"\t."}' \
		${clone}_${target}_${type}_${control}_control_SEACR.relaxed_rmBL.bed > ${clone}_${target}_${type}_${control}_control_SEACR.relaxed_rmBL.SAF
	"""

}

/*
 Step 15. Process TF bam files and control data to extract nucleosome-free region fragments
*/

process bam_to_NFR{
	publishDir "${resDir}/processed_beds", mode: 'copy'
	label 'bam_fragments'

	input:
	set val(clone), val(target), file(bam), val(type) from bam_TF_control

	when:
	params.dataTF == "yes"

	output:
	set val(clone), val(target), val(type), file("*.bed") into bed_NFR
	set val(clone), val(target), stdout into size_factors_NFR

	"""
	samtools sort -@ ${params.numCPUs} -n ${bam} -o ${bam.baseName}Name.bam
	bedtools bamtobed -bedpe -i ${bam.baseName}Name.bam | awk '\$1==\$4 && \$6-\$2 <= 120 {print\$1"\t"\$2"\t"\$6}' | \
		sort -k1,1 -k2,2n -k3,3n > ${bam.baseName}_NFR.bed
	echo 100000000/\$(wc -l ${bam.baseName}_NFR.bed | cut -d " " -f 1) | bc -l
	"""

}

/*
 Step 16. Generate normalized signal tracks using filtered pairs
*/

// Format channel with nucleosomal size factors to combine it with the bam files

size_factors_nucleosome
	.map{ row -> tuple(row[0] + "_" + row[1], row[2]) }
	.set { size_factors }

// Generate tracks with nucleosomal signal

process nucleosome_tracks{
	publishDir "${resDir}/signal/normalized_per_library/nucleosomal", mode: 'copy'
	label 'signal_tracks'

	input:
	set val(name), file(bam_cov), val(sizeFactor) from bam_coverage.combine(size_factors, by: 0)

	output:
	file("*.bw") into bigWig_nucleosome

	"""
	samtools index ${bam_cov}
	bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_nucleosomal.bw --outFileFormat bigwig --binSize 1 --normalizeUsing CPM \
		--effectiveGenomeSize 2407883318 --numberOfProcessors ${params.numCPUs_Dtools} --minFragmentLength 121 --maxFragmentLength 500 \
		--ignoreForNormalization chr1 chr8 chr16 chrX -e --scaleFactor ${sizeFactor}
	"""

}

// Generate channel with total and allelic data of TF and control libraries

G1_coverage
	.mix(G2_coverage)
	.map{ row -> tuple(row[0].split("_")[0] + "_" + row[0].split("_")[1], row[0].split("_")[2], row[1]) }
	.join(dataType2, by: [0,1])
	.filter{ it[3] == 'control' || 'TF' }
	.set { bam_allelic }

bam_allelic
	.mix(bam_NFR_signal)
	.set { bam_NFR }

// Generate tracks with NFR signal

process NFR_tracks{
	publishDir "${resDir}/signal/normalized_per_library/NFR", mode: 'copy'
	label 'signal_tracks'

	input:
	set val(clone), val(target), file(bam_cov), val(type), val(sizeFactor) from bam_NFR.combine(size_factors_NFR, by: [0,1])

	when:
	params.dataTF == "yes"
	
	output:
	file("*.bw") into bigWig_NFR

	"""
	samtools index ${bam_cov}
	bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_NFR.bw --outFileFormat bigwig --binSize 1 --normalizeUsing CPM \
		--effectiveGenomeSize 2407883318 --numberOfProcessors ${params.numCPUs_Dtools} --maxFragmentLength 120 \
		--ignoreForNormalization chr1 chr8 chr16 chrX -e --scaleFactor ${sizeFactor}
	"""

}

/*
 Step 17. Generate control normalized signal tracks of TFs (total signal only)
*/

process TF_norm{
	publishDir "${resDir}/signal/normalized_total_TF_signal", mode: 'copy'

	input:
	set val(clone), val(target), val(type), file(bed_target), val(control), file(bed_control) from bed_TF_tracks

	when:
	params.dataTF == "yes"

	output:
	file("*.bw") into bigWig_TF

	"""
	macs2 callpeak -f BEDPE -g mm --keep-dup all -B --nomodel --SPMR -n ${clone}_${target} -t ${bed_target} -c ${bed_control}
	macs2 bdgcmp -m FE -o ${clone}_${target}_${type}_${control}_normalized_FE.bdg -t ${clone}_${target}_treat_pileup.bdg \
		-c ${clone}_${target}_control_lambda.bdg
	sort -k1,1 -k2,2n ${clone}_${target}_${type}_${control}_normalized_FE.bdg > ${clone}_${target}_${type}_${control}_normalized_FE_sorted.bedgraph
	bedGraphToBigWig ${clone}_${target}_${type}_${control}_normalized_FE_sorted.bedgraph ${chrFile} \
		${clone}_${target}_${type}_${control}_normalized_FE.bw
	"""

}

/*
 Step 18. Assign peaks to genes based on genomic distance
*/

process peak_Annot{
	conda '/home/user/conda-envs/RpeakAnnotation'
	publishDir "${resDir}/qc/14_peak_annotation", mode: 'copy'
	label 'scR'

	input:
	set val(clone), val(target), val(caller), file(peaks) from peaks_annot_macs.mix(peaks_annot_seacr)

	output:
	file("${clone}_${target}*") into plots_annot

	"""
	Rscript ${plotPeakAnno} ${peaks} ${clone}_${target}_${caller}
	"""

}

/*
 Step 19. Calculate fraction of reads in peaks
*/

// Generate channels with peak and fragment data

peaks_macs
	.join(bed_frip_macs, by: [0,1])
	.set { peaks_frip_macs }

peaks_histone_seacr
	.join(bed_histone, by: [0,1])
	.set { peaks_frip_seacr }

// Calculate values

process peak_signal{
	publishDir "${resDir}/qc/15_fraction_reads_peaks", mode: 'copy'

	input:
	set val(clone), val(target), val(caller), file(peaks), file(fragments) from peaks_frip_macs.mix(peaks_frip_seacr)

	output:
	file('*FRiP.txt') into FRiP

	"""
	wc -l ${fragments} > ${clone}_${target}_total_fragments.txt
	bedtools sort -i ${peaks} | bedtools merge -i stdin | bedtools intersect -nonamecheck -u -a ${fragments} -b stdin | \
			wc -l > ${clone}_${target}_fragments_peaks.txt
	paste ${clone}_${target}_total_fragments.txt ${clone}_${target}_fragments_peaks.txt | \
			awk 'BEGIN{print"Total fragments\tFragments in peaks\tFRiP"};{print \$1"\t"\$3"\t"\$3/\$1}' > ${clone}_${target}_${caller}_FRiP.txt
	"""

}

/*
 Step 20. Read counting on peaks
*/

// Generate channel with peak information and bam files

bam_read_counts
	.map{ row -> tuple(row[0].split("_")[0] + "_" + row[0].split("_")[1], row[0].split("_")[2], row[1]) }
	.combine(dataType3, by: [0,1])
	.branch {
		control: it[3] == 'control'
		target: it[3] == 'TF' || 'histone'
	}
	.set { bam_annotated }

peaks_counting_macs.mix(peaks_counting_seacr).into { peaks_counting1; peaks_counting2 }

bam_annotated.target
	.combine(peaks_counting1, by: [0,1])
	.map{ row -> tuple(row[0], row[1], row[2], row[4], row[5], row[6]) }
	.set { peaks_counting_target }

bam_annotated.control
	.map{ row -> tuple(row[0], row[2], row[1]) }
	.combine(peaks_counting2, by: [0,2])
	.map{ row -> tuple(row[0], row[3], row[2], row[1], row[4], row[5]) }
	.mix(peaks_counting_target)
	.set { peaks_counting }

// Count number of fragments by peak

process read_counts{
	publishDir "${resDir}/quant/macs2_peaks", mode: 'copy', pattern: '*macs2.counts'
	publishDir "${resDir}/quant/SEACR_peaks", mode: 'copy', pattern: '*SEACR.counts'
	publishDir "${resDir}/qc/16_read_counting", mode: 'copy', pattern: '*.summary'

	input:
	set val(clone), val(target), file(bam), val(control), val(caller), file(peaks) from peaks_counting

	output:
	file("*macs2.counts") optional true into count_files_macs2
	file("*SEACR.counts") optional true into count_files_seacr
	file("*.summary") into counting_summary

	"""
	featureCounts -C -p -P -d 121 -D 500 -B -F SAF -a ${peaks} -T ${params.numCPUs} \
		-o ${bam.baseName}_countsOn_${clone}_${target}_${control}_control_${caller}.counts ${bam}
	"""

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Main results are saved in ${resDir}\n" : "There was an error during the execution, check log files." )
}


