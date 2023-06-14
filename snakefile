rule all:
    input:
       'trimSRR8234111.fastq.gz',
       'qc_output/trimSRR8234111_fastqc.html',
       'Tair10/trimSRR8234111.fastq.sam.gz',
       'Tair10/trimSRR8234111.fastq.bam',
       'Tair10/sortedtrimSRR8234111.fastq.bam',
       'Tair10/sortedtrimSRR8234111.fastq.bam.bai',
       'Rstudio/Scripts/Analyse_DEG.R',
       'Rstudio/DiffExprs/DEG_AllContrasts_Gitte_Dirk.csv',
       'atRTD3_29122021_index',
       'RTD3/SRR8234111/abundance.tsv'


rule trim:
    input:
        'SRR8234111.fastq.gz'
    output:
        'trimSRR8234111.fastq.gz'
    conda:
        'snakemake_envs/fastx_toolkit.yml'
    shell:
        'bash sript.sh'

rule quality_control:
    input:
        'trimSRR8234111.fastq'
    output:
        'qc_output/trimSRR8234111_fastqc.html'
    conda:
        'snakemake_envs/fastqc.yml'
    shell:
        'fastqc trim*.fastq -o qc_output/'
 

rule align_Tair10:
    input:
        'trimSRR8234111.fastq.gz'
    output:
        'Tair10/trimSRR8234111.fastq.sam.gz'
    conda:
        'snakemake_envs/hisat2.yml'
    shell:
        'bash allign.sh'
        

rule samtools_Tair10_bam:
    input:
        'Tair10/trimSRR8234111.fastq.sam'
    output:
        'Tair10/trimSRR8234111.fastq.bam'
    conda:
        'snakemake_envs/samtools.yml'
    shell:
        'bash Tair10/convert.sh'
     
rule sort_Tair10:
    input:
        'Tair10/trimSRR8234111.fastq.bam'
    output:
        'Tair10/sortedtrimSRR8234111.fastq.bam'
    conda:
        'snakemake_envs/samtools.yml'
    shell:
        'bash Tair10/sort.sh'

rule index_Tair10:
    input:
        'Tair10/sortedtrimSRR8234111.fastq.bam'
    output:
        'Tair10/sortedtrimSRR8234111.fastq.bam.bai'
    conda:
        'snakemake_envs/samtools.yml'
    shell:
        'bash Tair10/index.sh'

rule DEG:
    input:
        counts="Rstudio/Counts/Counts.txt",
        metadata="Rstudio/meta/Metadata_3_new.csv"
    output:
        deg_results="Rstudio/DiffExprs/DEG_AllContrasts_Gitte_Dirk.csv"
    shell:
        "Rscript {Rstudio/Scripts/functions.R} && Rscript {Rstudio/Scripts/Analyse_DEG.R} --counts {input.counts} --metadata {input.metadata} --output {output.deg_results}"

rule heatmap:
    input:
        deg_results="Rstudio/DiffExprs/DEG_AllContrasts_Gitte_Dirk.csv",
        clustering_script="Rstudio/Scripts/Clustering_Heatmap.R"
    output:
        heatmap_plot="Rstudio/heatmap/dif_SL_root_07071910_pam_heatmap.pdf"
    shell:
        "Rscript {input.clustering_script} --deg_results {input.deg_results} --output {output.heatmap_plot}"

rule kallisto_index_RTD3:
    input:
        'atRTD3_29122021.fa'
    output:
        'atRTD3_29122021_index'
    conda:
        'snakemake_envs/kallisto.yml'
    shell:
        'index --index=atRTD3_29122021_index atRTD3_29122021.fa'

rule kallisto_RTD3:
    input: 
        'atRTD3_29122021_index'
    output:
        'RTD3/SRR8234111/abundance.tsv'
    conda:
        'snakemake_envs/suppa.yml'
    shell:
        'RTD3/quant -i ../atRTD3_29122021_index -o SRR8234111 --single -b 100 -l 180 -s 20 -g atRTD3_TS_21Feb22_transfix.gtf ../trimSRR8234111.fastq.gz -t --pseudobam'
    
