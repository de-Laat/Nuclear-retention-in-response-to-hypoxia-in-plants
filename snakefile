rule all:
    input:
       'trimSRR8234111.fastq.gz',
       'qc_output/trimSRR8234111_fastqc.html',
       'Tair10/trimSRR8234111.fastq.sam.gz',
       'Tair10/trimSRR8234111.fastq.bam',
       'Tair10/sortedtrimSRR8234111.fastq.bam',
       'Tair10/sortedtrimSRR8234111.fastq.bam.bai',
       "Rstudio/Counts/Counts.txt",
       "Rstudio/meta/Metadata_3_new.csv",
       "Rstudio/DiffExprs/DEG_AllContrasts_Gitte_Dirk.csv",
       "Rstudio/dif_SL_root_070719_k10_GO_IN.txt",
       "Rstudio/heatmap/dif_SL_root_07071910_pam_heatmap.pdf",
       "Rstudio/GOslimbarplotBP_1.pdf",
       'atRTD3_29122021_index',
       'RTD3/SRR8234111/abundance.tsv',
       'suppa_ioe/suppa_ioe_RI_strict.ioe',
       'test_suppa/all_reps_con1_tpms.tpm',
       'splice_events/nuc_2ns_RI.psi'
       

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
        deg_results1="Rstudio/DiffExprs/DEG_AllContrasts_Gitte_Dirk.csv",
        deg_results2="Rstudio/dif_SL_root_070719_k10_GO_IN.txt"
    shell:
        "Rscript {Rstudio/Scripts/functions.R} && Rscript {Rstudio/Scripts/Analyse_DEG.R} --counts {input.counts} --metadata {input.metadata} --output {output.deg_results}"

rule heatmap:
    input:
        deg_results1="Rstudio/DiffExprs/DEG_AllContrasts_Gitte_Dirk.csv"
    output:
        heatmap_plot="Rstudio/heatmap/dif_SL_root_07071910_pam_heatmap.pdf"
    shell:
        "Rscript {Rstudio/Scripts/Clustering_Heatmap.R} --deg_results1 {input.deg_results1} --output {output.heatmap_plot}"

rule GEA:
    input:
        deg_results2="Rstudio/dif_SL_root_070719_k10_GO_IN.txt"
    output:
        barplot="Rstudio/GOslimbarplotBP_1.pdf"
    shell:
        "Rscript {Rstudio/Scripts/GEA.R} --deg_results2 {input.deg_results2} --output {output.barplot}

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
        'snakemake_envs/kallisto.yml'
    shell:
        'RTD3/quant -i ../atRTD3_29122021_index -o SRR8234111 --single -b 100 -l 180 -s 20 -g atRTD3_TS_21Feb22_transfix.gtf ../trimSRR8234111.fastq.gz -t --pseudobam'
    
rule suppa_generate_events:
    input:
        'atRTD3_29122021_index'
    output:
        'suppa_ioe/suppa_ioe_RI_strict.ioe'
    conda:
        'snakemake_envs/suppa.yml'
    shell:
        'suppa.py generateEvents -i ../atRTD3_TS_21Feb22_transfix.gtf -o suppa_ioe -f ioe -e SE SS MX RI FL'

rule suppa_join_files:
    input:
        'RTD3/SRR8234111/abundance.tsv',
        'RTD3/SRR8234115/abundance.tsv',
        'RTD3/SRR8234119/abundance.tsv'
    output:
        'test_suppa/all_reps_con1_tpms.tpm'
    conda:
        'snakemake_envs/suppa.yml'
    shell:
        'suppa.py joinFiles -f tpm -i abundanceSRR8234111.txt abundanceSRR8234115.txt abundanceSRR8234119.txt -o all_reps_con1_tpms '

rule suppa_psiPerEvent:
    input:
        'test_suppa/all_reps_con1_tpms.tpm',
        'suppa_ioe/suppa_ioe_RI_strict.ioe'
    output:
        'splice_events/nuc_2ns_RI.psi'
    conda:
        'snakemake_envs/suppa.yml'
    shell:
        'suppa.py psiPerEvent -e reps_nuc_2ns_tpm/all_reps_con1_tpms.tpm -i suppa_ioe/suppa_ioe_RI_strict.ioe -o splice_events/nuc_2ns_RI'

rule diff_suppa:
    input:
        'splice_events/nuc_2ns_RI.psi',
        'splice_events/nuc_9ns/nuc_9ns_RI.psi',
        'suppa_ioe/suppa_ioe_RI_strict.ioe',
        'test_suppa/all_reps_con1_tpms.tpm',
        'reps_nuc_2hs_tpm/all_reps_con2_tpms.tpm'
    output:
        'diffsplice_events2/2_nuc_RI.dpsi.temp.0'
    conda:
        'snakemake_envs/suppa.yml'
    shell:
        'suppa.py diffSplice -m empirical -i suppa_ioe/suppa_ioe_RI_strict.ioe -p splice_events/nuc_9hs/nuc_9hs_RI.psi splice_events/nuc_9ns/nuc_9ns_RI.psi -e reps_nuc_9hs_tpm/all_reps_con4_tpms.tpm reps_nuc_9ns_tpm/all_reps_con3_tpms.tpm --save_tpm_events -l 0.05 -pa -gc -c -o diffsplice_events2/9_NUC_RI'
    
