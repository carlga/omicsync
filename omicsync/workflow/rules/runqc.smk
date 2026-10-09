"""
================================================================================
RUNQC MODULE - Quality Control and Preprocessing Workflow
================================================================================
Performs comprehensive QC on sequencing data including:
- Subsampling and quality assessment (FastQC)
- Read preprocessing (fastp)
- Optional biotype composition analysis (RSeQC, HISAT2, featureCounts)
- Optional species composition analysis (Kraken2, Bracken, Krona)
- MultiQC report generation
================================================================================
"""

import os
import re
import glob
from pathlib import Path


# Setup workflow environment
workflow_dir = Path(workflow.snakefile).parent.absolute()
sys.path.insert(0, str(workflow_dir.parent))


# Import utilities
from utils import check_strandedness, plot_biotype_composition
from utils.workflow_utils import get_timestamp


# Module identifier
MODULE = 'runqc'
TIMESTAMP = get_timestamp()



# ================================================================================
# CONFIGURATION AND SETUP
# ================================================================================

# Set variables (quiet parsing - details shown in onstart)
PROJID = os.path.basename(config['runqc']['indir'])
FILES = glob.glob(os.path.join(config['runqc']['indir'], '*.f*q.gz'))
EXT = re.search(r'[-\w]+(.f(ast)?q.gz)', FILES[0]).group(1)

try:
    NAMES = list(set([re.search(r'([-\w]+)(?=_[12])_[12].f(ast)?q.gz', x).group(1) for x in FILES]))
except AttributeError:
    NAMES = list(set([re.search(r'([-\w]+).f(ast)?q.gz', x).group(1) for x in FILES]))

PAIRED = False if len(NAMES) == len(FILES) else True
READS = ['_1','_2'] if PAIRED else ['']

# Set cores and opt
AVAIL_THREADS = config['runqc']['jobs']
try:
    FASTP_THREADS = int(re.search(r'(-w|--thread) (\d+)', config['runqc']['fastp_opt']).group(2))
    config['runqc']['fastp_opt'] = re.sub(r'(-w|--thread) (\d+) ', '', config['runqc']['fastp_opt'])
except AttributeError:
    FASTP_THREADS = 3
FASTP_THREADS = AVAIL_THREADS if FASTP_THREADS > AVAIL_THREADS else FASTP_THREADS

# Check for optional analyses
if config["runqc"]["hisat2_ref"] and config["runqc"]["gtf"]:
    HISAT2_THREADS = FASTP_THREADS

if config["runqc"]["kraken2_db"]:
    KRAKEN2_THREADS = AVAIL_THREADS
    BRACKEN_DB = glob.glob(os.path.join(config['runqc']['kraken2_db'], '*mers.kmer_distrib'))
    if len(BRACKEN_DB) == 1:
        BRACKEN_READ_LEN = int(re.search(r'database(\d+)mers.kmer_distrib', BRACKEN_DB[0]).group(1))
    else:
        BRACKEN_READ_LEN = 100



# ================================================================================
# WORKFLOW INITIALIZATION
# ================================================================================

onstart:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Quality Control Workflow'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Project: {PROJID}'
        echo '[{MODULE}] Samples: {len(NAMES)}'
        echo '[{MODULE}] Layout: {'Paired-end' if PAIRED else 'Single-end'}'
        echo '[{MODULE}] Sampling: {'Full dataset' if config['runqc']['number'] == 1 else f"{config['runqc']['number']}{'%' if config['runqc']['number'] < 1 else ' reads'}"}'
        echo ''
        echo '[{MODULE}] Enabled analyses:'
        echo '[{MODULE}]   • FastQC (raw reads)'
        echo '[{MODULE}]   • Read preprocessing (fastp)'
        echo '[{MODULE}]   • FastQC (preprocessed reads)'
        if [ -n "{config['runqc'].get('hisat2_ref', '')}" ] && [ -n "{config['runqc'].get('gtf', '')}" ]; then
            echo '[{MODULE}]   • Biotype composition (RSeQC/HISAT2/featureCounts)'
        fi
        if [ -n "{config['runqc'].get('kraken2_db', '')}" ]; then
            echo '[{MODULE}]   • Species composition (Kraken2/Bracken/Krona)'
        fi
        echo '[{MODULE}]   • MultiQC report'
        echo ''
        echo '[{MODULE}] Threads: {AVAIL_THREADS}'
        echo '[{MODULE}] Timestamp: {TIMESTAMP}'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo ''
    """)



# ================================================================================
# TARGET RULES
# ================================================================================

# Set targets
RUN_RULES = expand(
    ["fq_raw_sampled/{name}_subs{read}.fastq.gz",
    "fastqc_raw/{name}_subs{read}_fastqc.html",
    "fq_preprocessed/{name}_fastp.json",
    "fastqc_preprocessed/{name}_subs_preproc{read}_fastqc.html"],
    name=NAMES, read=READS
)

RUNQC_BIOTYPE_COMPOSITION = RUNQC_SPECIES_COMPOSITION = None

if config["runqc"]["hisat2_ref"] and config["runqc"]["gtf"]:
    GTF_BASENAME = os.path.splitext(os.path.basename(config["runqc"]["gtf"]))[0]
    RUNQC_BIOTYPE_COMPOSITION = expand(
        ["rseqc/{name}_rseqc_inferexperiment.txt",
        "rseqc/{name}_rseqc_read_distribution.txt",
        "hisat2_aligned/{name}.bam",
        "hisat2_aligned_featurecounts/" + PROJID + "_subs_featurecounts.txt"],
        name=NAMES, read=READS
    )
    RUN_RULES.append(RUNQC_BIOTYPE_COMPOSITION)

if config["runqc"]["kraken2_db"]:
    RUNQC_SPECIES_COMPOSITION = expand(
        ["kraken2/{name}_kraken2_report.txt",
        "bracken/{name}_kraken2_bracken_report.txt",
        "krona/" + PROJID + "_species_composition.html"],
        name=NAMES, read=READS
    )
    RUN_RULES.append(RUNQC_SPECIES_COMPOSITION)



rule all:
    input:
        RUN_RULES + ["multiqc/multiqc_report.html"]



# ================================================================================
# FASTQ SAMPLING AND QUALITY CONTROL
# ================================================================================


rule subset_fq:
    """
    Subsamples raw fastq files to desired read number/proportion.
    Supports three modes:
      - number < 1: sample proportion (e.g., 0.1 = 10%)
      - number = 1: copy entire file
      - number > 1: sample exact number of reads
    """
    input:
        fq = config['runqc']['indir'] + "/{name}{read}" + EXT
    output:
        fq_subs = "fq_raw_sampled/{name}_subs{read,.*}.fastq.gz"
    params:
        number = config['runqc']['number'],
        seed = config['runqc']['seed']
    log: "logs/" + TIMESTAMP + "_{name}{read}_subset_fq.log"
    threads: 1
    run:
        shell(f"echo '[{MODULE}] Sampling {{wildcards.name}}{{wildcards.read}}{EXT} (p={{params.number}})...' > {{log}}")

        if params['number'] < 1:
            shell(
                f'''
                seqkit sample \
                    -j {{threads}} \
                    -p {{params.number}} \
                    -s {{params.seed}} \
                    -o {{output.fq_subs}} \
                    {{input.fq}} \
                    &>> {{log}}
                '''
            )
        elif params['number'] == 1:
            shell("cp {input.fq} {output.fq_subs}")
        elif params['number'] > 1:
            shell(
                '''
                seqkit sample \
                    -j {threads} \
                    -p 1 \
                    -s {params.seed} \
                    -o {output.fq_subs}.tmp \
                    {input.fq} \
                    2>> {log}
                
                seqkit head \
                    -n {params.number} \
                    -o {output.fq_subs} \
                    {output.fq_subs}.tmp \
                    &>> {log}

                rm {output.fq_subs}.tmp
                '''
            )

        shell(f"echo '[{MODULE}] ✓ Sampling complete for {{wildcards.name}}{{wildcards.read}}{EXT}' >> {{log}}")



rule fastqc_raw:
    """
    Runs FastQC on sampled raw fastq files.
    """
    input:
        fq_files = expand(
            "fq_raw_sampled/{name}_subs{read}.fastq.gz", 
            name=NAMES, read=READS
        )
    output:
        qc_reports = expand(
            "fastqc_raw/{name}_subs{read}_fastqc.html", 
            name=NAMES, read=READS
        )
    params:
        outdir = "fastqc_raw"
    log: "logs/" + TIMESTAMP + "_fastqc_raw.log"
    threads: AVAIL_THREADS
    shell:
        f'''
        echo '[{MODULE}] Running FastQC on raw fastq files...' > {{log}}

        fastqc \
            --threads {{threads}} \
            --outdir {{params.outdir}} \
            {{input.fq_files}} \
            &>> {{log}}

        echo '[{MODULE}] ✓ FastQC complete for raw files' >> {{log}}
        '''



# ================================================================================
# READ PREPROCESSING
# ================================================================================


if PAIRED:
    rule preprocess_fq:
        """
        Performs preprocessing on paired-end raw fastq files using fastp.
        Includes adapter trimming, quality filtering, and length filtering.
        """
        input:
            fq1 = "fq_raw_sampled/{name}_subs_1.fastq.gz",
            fq2 = "fq_raw_sampled/{name}_subs_2.fastq.gz"
        output:
            fq1_preproc = "fq_preprocessed/{name}_subs_preproc_1.fastq.gz",
            fq2_preproc = "fq_preprocessed/{name}_subs_preproc_2.fastq.gz",
            json = "fq_preprocessed/{name}_fastp.json",
            html = "fq_preprocessed/{name}_fastp.html"
        params:
            fastp_opt = config['runqc']['fastp_opt']
        log: "logs/" + TIMESTAMP + "_{name}_preprocess_fq.log"
        threads: FASTP_THREADS
        shell:
            f'''
            echo '[{MODULE}] Processing paired-end reads from {{wildcards.name}}...' > {{log}}

            fastp \
                --thread {{threads}} \
                {{params.fastp_opt}} \
                -i {{input.fq1}} \
                -o {{output.fq1_preproc}} \
                -I {{input.fq2}} \
                -O {{output.fq2_preproc}} \
                -j {{output.json}} \
                -h {{output.html}} \
                &>> {{log}}

            echo '[{MODULE}] ✓ Preprocessing complete for {{wildcards.name}}' >> {{log}}
            '''
else:
    rule preprocess_fq:
        """
        Performs preprocessing on single-end raw fastq files using fastp.
        Includes adapter trimming, quality filtering, and length filtering.
        """
        input:
            fq = "fq_raw_sampled/{name}_subs.fastq.gz"
        output:
            fq_preproc = "fq_preprocessed/{name}_subs_preproc.fastq.gz",
            json = "fq_preprocessed/{name}_fastp.json",
            html = "fq_preprocessed/{name}_fastp.html"
        params:
            fastp_opt = config['runqc']['fastp_opt']
        log: "logs/" + TIMESTAMP + "_{name}_preprocess_fq.log"
        threads: FASTP_THREADS
        shell:
            f'''
            echo '[{MODULE}] Processing single-end reads from {{wildcards.name}}...' > {{log}}

            fastp \
                --thread {{threads}} \
                {{params.fastp_opt}} \
                -i {{input.fq}} \
                -o {{output.fq_preproc}} \
                -j {{output.json}} \
                -h {{output.html}} \
                &>> {{log}}

            echo '[{MODULE}] ✓ Preprocessing complete for {{wildcards.name}}' >> {{log}}
            '''



rule fastqc_preprocessed:
    """
    Runs FastQC on sampled and preprocessed fastq files.
    """
    input:
        fq_files = expand(
            "fq_preprocessed/{name}_subs_preproc{read}.fastq.gz",
            name=NAMES, read=READS
        )
    output:
        qc_reports = expand(
            "fastqc_preprocessed/{name}_subs_preproc{read}_fastqc.html",
            name=NAMES, read=READS
        )
    params:
        outdir="fastqc_preprocessed"
    log: "logs/" + TIMESTAMP + "_fastqc_preprocessed.log"
    threads: AVAIL_THREADS
    shell:
        f'''
        echo '[{MODULE}] Running FastQC on preprocessed fastq files...' > {{log}}

        fastqc \
            --threads {{threads}} \
            --outdir {{params.outdir}} \
            {{input.fq_files}} \
            &>> {{log}}

        echo '[{MODULE}] ✓ FastQC complete for preprocessed files' >> {{log}}
        '''



# ================================================================================
# BIOTYPE COMPOSITION ANALYSIS (OPTIONAL)
# ================================================================================


if RUNQC_BIOTYPE_COMPOSITION:

    rule rseqc_bed:
        """
        Converts GTF annotation to BED format for RSeQC analysis.
        Uses UCSC tools gtfToGenePred and genePredToBed.
        """
        input:
            gtf = config["runqc"]["gtf"]
        output:
            bed = "rseqc/" + GTF_BASENAME + "_transcripts.bed"
        params:
            outdir = "rseqc"
        log: "logs/" + TIMESTAMP + "_rseqc_bed.log"
        threads: 1
        shell:
            f'''
            echo '[{MODULE}] Converting GTF to BED for RSeQC analysis...' > {{log}}

            gtfToGenePred \
                {{input.gtf}} \
                {{params.outdir}}/transcript_genepred.tmp \
                &>> {{log}}
            
            genePredToBed \
                {{params.outdir}}/transcript_genepred.tmp \
                {{output.bed}} \
                &>> {{log}}

            rm {{params.outdir}}/transcript_genepred.tmp

            echo '[{MODULE}] ✓ GTF to BED conversion complete' >> {{log}}
            '''
    
    if PAIRED:
    
        rule rseqc_infer_experiment:
            """
            Infers strandedness of paired-end sequencing experiment using RSeQC.
            Performs quick alignment with HISAT2 to determine library strandedness.
            """
            input:
                fq1 = "fq_preprocessed/{name}_subs_preproc_1.fastq.gz",
                fq2 = "fq_preprocessed/{name}_subs_preproc_2.fastq.gz",
                bed = "rseqc/" + GTF_BASENAME + "_transcripts.bed"
            output:
                infer_experiment = "rseqc/{name}_rseqc_inferexperiment.txt"
            params:
                hisat2_ref = config['runqc']['hisat2_ref'],
                outdir = "rseqc"
            log: "logs/" + TIMESTAMP + "_{name}_rseqc_infer_experiment.log"
            threads: HISAT2_THREADS
            shell:
                f'''
                echo '[{MODULE}] Inferring strandedness for {{wildcards.name}}...' > {{log}}

                # Align reads to reference genome
                hisat2 \
                    --threads {{threads}} \
                    -1 {{input.fq1}} \
                    -2 {{input.fq2}} \
                    -x {{params.hisat2_ref}} \
                    -S {{params.outdir}}/{{wildcards.name}}.sam \
                    &>> {{log}}

                # Convert to BAM, filter unmapped reads (0x4, 0x8, 0x100)
                samtools view \
                    -@ {{threads}} \
                    -bS \
                    -F 268 \
                    {{params.outdir}}/{{wildcards.name}}.sam \
                    > {{params.outdir}}/{{wildcards.name}}.bam \
                    2>> {{log}}
                
                samtools sort \
                    -@ {{threads}} \
                    {{params.outdir}}/{{wildcards.name}}.bam \
                    -o {{params.outdir}}/{{wildcards.name}}.bam \
                    &>> {{log}}
                
                samtools index \
                    -@ {{threads}} \
                    {{params.outdir}}/{{wildcards.name}}.bam \
                    &>> {{log}}

                # Infer experiment strandedness
                infer_experiment.py \
                    -i {{params.outdir}}/{{wildcards.name}}.bam \
                    -r {{input.bed}} \
                    > {{output.infer_experiment}} \
                    2>> {{log}}

                # Cleanup temporary files
                rm {{params.outdir}}/{{wildcards.name}}.sam
                rm {{params.outdir}}/{{wildcards.name}}.bam
                rm {{params.outdir}}/{{wildcards.name}}.bam.bai

                echo '[{MODULE}] ✓ Strandedness inference complete for {{wildcards.name}}' >> {{log}}
                '''
        
        rule hisat2_align:
            """
            Aligns paired-end reads to reference genome using HISAT2.
            Uses inferred strandedness for strand-specific alignment.
            """
            input:
                fq1 = "fq_preprocessed/{name}_subs_preproc_1.fastq.gz",
                fq2 = "fq_preprocessed/{name}_subs_preproc_2.fastq.gz",
                infer_experiment = "rseqc/{name}_rseqc_inferexperiment.txt"
            output:
                sam = temp("hisat2_aligned/{name}.sam"),
                bam = "hisat2_aligned/{name}.bam",
                bam_idx = "hisat2_aligned/{name}.bam.bai",
                summary = "hisat2_aligned/{name}_hisat2_summary.txt",
                stats = "hisat2_aligned/{name}_stats.txt"
            params:
                hisat2_ref = config['runqc']['hisat2_ref'],
                outdir = "hisat2_aligned"
            log: "logs/" + TIMESTAMP + "_{name}_hisat2_align.log"
            threads: HISAT2_THREADS
            run:
                strandedness = check_strandedness(input['infer_experiment'], PAIRED)

                shell(f"echo '[{MODULE}] Inferred strandedness for {{wildcards.name}}: {strandedness}' > {{log}}")

                hisat2_strandedness = ''
                if strandedness == 'fwd':
                    hisat2_strandedness = '--rna-strandness FR'
                elif strandedness == 'rev':
                    hisat2_strandedness = '--rna-strandness RF'
                
                shell(
                    f'''
                    echo '[{MODULE}] Aligning {{wildcards.name}} with HISAT2...' >> {{log}}

                    hisat2 \
                        --threads {{threads}} \
                        {hisat2_strandedness} \
                        -1 {{input.fq1}} \
                        -2 {{input.fq2}} \
                        -x {{params.hisat2_ref}} \
                        --summary-file {{output.summary}} \
                        -S {{output.sam}} \
                        &>> {{log}}

                    echo '[{MODULE}] Processing alignment results...' >> {{log}}

                    # Convert to BAM, filter unmapped reads (0x4, 0x8, 0x100)
                    samtools view \
                        -@ {{threads}} \
                        -bS \
                        -F 268 \
                        {{output.sam}} \
                        > {{output.bam}} \
                        2>> {{log}}
                    
                    samtools sort \
                        -@ {{threads}} \
                        {{output.bam}} \
                        -o {{output.bam}} \
                        &>> {{log}}
                    
                    samtools index \
                        -@ {{threads}} \
                        {{output.bam}} \
                        &>> {{log}}

                    samtools flagstat \
                        -@ {{threads}} \
                        {{output.bam}} \
                        > {{output.stats}} \
                        2>> {{log}}

                    echo '[{MODULE}] ✓ Alignment complete for {{wildcards.name}}' >> {{log}}
                    '''
                )
    
    else:

        rule rseqc_infer_experiment:
            """
            Infers strandedness of single-end sequencing experiment using RSeQC.
            Performs quick alignment with HISAT2 to determine library strandedness.
            """
            input:
                fq = "fq_preprocessed/{name}_subs_preproc.fastq.gz",
                bed = "rseqc/" + GTF_BASENAME + "_transcripts.bed"
            output:
                infer_experiment = "rseqc/{name}_rseqc_inferexperiment.txt"
            params:
                hisat2_ref = config['runqc']['hisat2_ref'],
                outdir = "rseqc"
            log: "logs/" + TIMESTAMP + "_{name}_rseqc_infer_experiment.log"
            threads: HISAT2_THREADS
            shell:
                f'''
                echo '[{MODULE}] Inferring strandedness for {{wildcards.name}}...' > {{log}}

                # Align reads to reference genome
                hisat2 \
                    --threads {{threads}} \
                    -U {{input.fq}} \
                    -x {{params.hisat2_ref}} \
                    -S {{params.outdir}}/{{wildcards.name}}.sam \
                    &>> {{log}}

                # Convert to BAM, filter unmapped reads (0x4, 0x100)
                samtools view \
                    -@ {{threads}} \
                    -bS \
                    -F 260 \
                    {{params.outdir}}/{{wildcards.name}}.sam \
                    > {{params.outdir}}/{{wildcards.name}}.bam \
                    2>> {{log}}
                
                samtools sort \
                    -@ {{threads}} \
                    {{params.outdir}}/{{wildcards.name}}.bam \
                    -o {{params.outdir}}/{{wildcards.name}}.bam \
                    &>> {{log}}
                
                samtools index \
                    -@ {{threads}} \
                    {{params.outdir}}/{{wildcards.name}}.bam \
                    &>> {{log}}

                # Infer experiment strandedness
                infer_experiment.py \
                    -i {{params.outdir}}/{{wildcards.name}}.bam \
                    -r {{input.bed}} \
                    > {{output.infer_experiment}} \
                    2>> {{log}}

                # Cleanup temporary files
                rm {{params.outdir}}/{{wildcards.name}}.sam
                rm {{params.outdir}}/{{wildcards.name}}.bam
                rm {{params.outdir}}/{{wildcards.name}}.bam.bai

                echo '[{MODULE}] ✓ Strandedness inference complete for {{wildcards.name}}' >> {{log}}
                '''
        
        rule hisat2_align:
            """
            Aligns single-end reads to reference genome using HISAT2.
            Uses inferred strandedness for strand-specific alignment.
            """
            input:
                fq = "fq_preprocessed/{name}_subs_preproc.fastq.gz",
                infer_experiment = "rseqc/{name}_rseqc_inferexperiment.txt"
            output:
                sam = temp("hisat2_aligned/{name}.sam"),
                bam = "hisat2_aligned/{name}.bam",
                bam_idx = "hisat2_aligned/{name}.bam.bai",
                summary = "hisat2_aligned/{name}_hisat2_summary.txt",
                stats = "hisat2_aligned/{name}_stats.txt"
            params:
                hisat2_ref = config['runqc']['hisat2_ref'],
                outdir = "hisat2_aligned"
            log: "logs/" + TIMESTAMP + "_{name}_hisat2_align.log"
            threads: HISAT2_THREADS
            run:
                strandedness = check_strandedness(input['infer_experiment'], PAIRED)

                shell(f"echo '[{MODULE}] Inferred strandedness for {{wildcards.name}}: {strandedness}' > {{log}}")

                hisat2_strandedness = ''
                if strandedness == 'fwd':
                    hisat2_strandedness = '--rna-strandness F'
                elif strandedness == 'rev':
                    hisat2_strandedness = '--rna-strandness R'

                shell(
                    f'''
                    echo '[{MODULE}] Aligning {{wildcards.name}} with HISAT2...' >> {{log}}

                    hisat2 \
                        --threads {{threads}} \
                        {hisat2_strandedness} \
                        -U {{input.fq}} \
                        -x {{params.hisat2_ref}} \
                        --summary-file {{output.summary}} \
                        -S {{output.sam}} \
                        &>> {{log}}

                    echo '[{MODULE}] Processing alignment results...' >> {{log}}

                    # Convert to BAM, filter unmapped reads (0x4, 0x100)
                    samtools view \
                        -@ {{threads}} \
                        -bS \
                        -F 260 \
                        {{output.sam}} \
                        > {{output.bam}} \
                        2>> {{log}}
                    
                    samtools sort \
                        -@ {{threads}} \
                        {{output.bam}} \
                        -o {{output.bam}} \
                        &>> {{log}}
                    
                    samtools index \
                        -@ {{threads}} \
                        {{output.bam}} \
                        &>> {{log}}

                    samtools flagstat \
                        -@ {{threads}} \
                        {{output.bam}} \
                        > {{output.stats}} \
                        2>> {{log}}

                    echo '[{MODULE}] ✓ Alignment complete for {{wildcards.name}}' >> {{log}}
                    '''
                )
    
    rule rseqc_read_distribution:
        """
        Analyzes distribution of mapped reads across genomic features using RSeQC.
        Reports read counts in CDS, UTRs, introns, intergenic regions, etc.
        """
        input:
            bam = "hisat2_aligned/{name}.bam",
            bed = "rseqc/" + GTF_BASENAME + "_transcripts.bed"
        output:
            read_distribution = "rseqc/{name}_rseqc_read_distribution.txt"
        log: "logs/" + TIMESTAMP + "_{name}_rseqc_read_distribution.log"
        threads: 1
        shell:
            f'''
            echo '[{MODULE}] Measuring read distribution for {{wildcards.name}}...' > {{log}}

            read_distribution.py \
                -i {{input.bam}} \
                -r {{input.bed}} \
                > {{output.read_distribution}} \
                2>> {{log}}

            echo '[{MODULE}] ✓ Read distribution analysis complete for {{wildcards.name}}' >> {{log}}
            '''

    rule feature_counts:
        """
        Quantifies gene biotypes using featureCounts from Subread package.
        Generates counts at biotype level (protein_coding, lncRNA, etc.).
        """
        input:
            bams = expand(["hisat2_aligned/{name}.bam"], name=NAMES),
            gtf = config['runqc']['gtf'],
            infer_experiment = "rseqc/" + NAMES[0] + "_rseqc_inferexperiment.txt"
        output:
            featurecounts = "hisat2_aligned_featurecounts/" + PROJID + "_subs_featurecounts.txt",
            biotype_composition = "hisat2_aligned_featurecounts/" + PROJID + "_biotype_composition.html"
        params:
            biotype = "gene_biotype"
        log: "logs/" + TIMESTAMP + "_feature_counts.log"
        threads: AVAIL_THREADS
        run:
            pe = '-p --countReadPairs -C -B' if PAIRED else ''

            strandedness = check_strandedness(input['infer_experiment'], PAIRED)

            shell(f"echo '[{MODULE}] Using inferred strandedness for {NAMES[0]}: {strandedness}' > {{log}}")

            featurecounts_strandedness = '-s 0'
            if strandedness == 'fwd':
                featurecounts_strandedness = '-s 1'
            elif strandedness == 'rev':
                featurecounts_strandedness = '-s 2'

            shell(
                f'''
                echo '[{MODULE}] Quantifying biotypes with featureCounts...' >> {{log}}

                featureCounts \
                    -T {{threads}} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -t exon \
                    -g {{params.biotype}} \
                    -a {{input.gtf}} \
                    -o {{output.featurecounts}} \
                    {{input.bams}} \
                    &>> {{log}}

                echo '[{MODULE}] ✓ Biotype quantification complete' >> {{log}}
                '''
            )

            plot_biotype_composition(output.featurecounts, output.biotype_composition)



# ================================================================================
# SPECIES COMPOSITION ANALYSIS (OPTIONAL)
# ================================================================================


if RUNQC_SPECIES_COMPOSITION:

    if PAIRED:

        rule kraken2:
            """
            Performs taxonomic classification of paired-end reads using Kraken2.
            Identifies species composition in the sample.
            """
            input:
                fq1 = "fq_preprocessed/{name}_subs_preproc_1.fastq.gz",
                fq2 = "fq_preprocessed/{name}_subs_preproc_2.fastq.gz"
            output:
                out = "kraken2/{name}_kraken2.txt",
                report = "kraken2/{name}_kraken2_report.txt"
            params:
                db = config['runqc']['kraken2_db']
            log: "logs/" + TIMESTAMP + "_{name}_kraken2.log"
            threads: KRAKEN2_THREADS
            shell:
                f'''
                echo '[{MODULE}] Running Kraken2 for {{wildcards.name}}...' > {{log}}

                kraken2 \
                    --threads {{threads}} \
                    --paired \
                    --gzip-compressed \
                    --db {{params.db}} \
                    --report {{output.report}} \
                    --output {{output.out}} \
                    {{input.fq1}} \
                    {{input.fq2}} \
                    &>> {{log}}

                echo '[{MODULE}] ✓ Kraken2 classification complete for {{wildcards.name}}' >> {{log}}
                '''

    else:

        rule kraken2:
            """
            Performs taxonomic classification of single-end reads using Kraken2.
            Identifies species composition in the sample.
            """
            input:
                fq = "fq_preprocessed/{name}_subs_preproc.fastq.gz"
            output:
                out = "kraken2/{name}_kraken2.txt",
                report = "kraken2/{name}_kraken2_report.txt"
            params:
                db = config['runqc']['kraken2_db']
            log: "logs/" + TIMESTAMP + "_{name}_kraken2.log"
            threads: KRAKEN2_THREADS
            shell:
                f'''
                echo '[{MODULE}] Running Kraken2 for {{wildcards.name}}...' > {{log}}

                kraken2 \
                    --threads {{threads}} \
                    --gzip-compressed \
                    --db {{params.db}} \
                    --report {{output.report}} \
                    --output {{output.out}} \
                    {{input.fq}} \
                    &>> {{log}}

                echo '[{MODULE}] ✓ Kraken2 classification complete for {{wildcards.name}}' >> {{log}}
                '''
    
    rule bracken:
        """
        Estimates species-level abundance from Kraken2 results using Bracken.
        Provides more accurate abundance estimates at species level.
        """
        input:
            k2report = "kraken2/{name}_kraken2_report.txt"
        output:
            out = "bracken/{name}_kraken2_bracken.txt",
            report = "bracken/{name}_kraken2_bracken_report.txt"
        params:
            db = config["runqc"]["kraken2_db"],
            bracken_read_len = BRACKEN_READ_LEN
        log: f"logs/{TIMESTAMP}_{{name}}_bracken.log"
        threads: KRAKEN2_THREADS
        shell:
            f'''
            echo '[{MODULE}] Running Bracken on {{wildcards.name}}...' | tee {{log}}
            
            # Run Bracken, but don't fail if it doesn't work
            if bracken \
                -r {{params.bracken_read_len}} \
                -l S \
                -t 5 \
                -d {{params.db}} \
                -i {{input.k2report}} \
                -w {{output.report}} \
                -o {{output.out}} \
                2>> {{log}}; then
                echo '[{MODULE}] ✓ Bracken abundance estimation complete for {{wildcards.name}}' | tee -a {{log}}
            else
                echo '[{MODULE}] ⚠ Bracken failed - creating empty output files' | tee -a {{log}}
                echo "# Bracken failed - no results" > {{output.out}}
                echo "# Bracken failed - no results" > {{output.report}}
            fi
            '''

    rule krona:
        """
        Generates interactive Krona plots for species composition visualization.
        Creates hierarchical taxonomy plot from Bracken results.
        Skips visualization if Bracken failed for all samples.
        """
        input:
            reports = expand(["bracken/{name}_kraken2_bracken_report.txt"], name=NAMES)
        output:
            krona_html = "krona/" + PROJID + "_species_composition.html"
        params:
            outdir = "krona"
        log: "logs/" + TIMESTAMP + "_krona.log"
        threads: 1
        shell:
            f'''
            echo '[{MODULE}] Generating Krona plots...' > {{log}}
            
            # Check if any Bracken reports have actual results
            has_results=0
            for report in {{input.reports}}; do
                if grep -v "^#" "$report" | grep -q "[0-9]"; then
                    has_results=1
                    break
                fi
            done
            
            if [ $has_results -eq 0 ]; then
                echo '[{MODULE}] ⚠ All Bracken reports are empty - creating placeholder Krona output' | tee -a {{log}}
                echo "<html><body><h1>No taxonomic classification results</h1>" > {{output.krona_html}}
                echo "<p>Bracken failed for all samples. This may be due to:</p>" >> {{output.krona_html}}
                echo "<ul><li>Database mismatch (e.g., using 16S database on RNA-seq data)</li>" >> {{output.krona_html}}
                echo "<li>No classified reads from Kraken2</li>" >> {{output.krona_html}}
                echo "<li>Read length incompatible with Bracken database</li></ul>" >> {{output.krona_html}}
                echo "</body></html>" >> {{output.krona_html}}
                exit 0
            fi

            # Prepare Bracken reports for Krona (only process non-empty ones)
            printf "%s\\n" {" ".join(NAMES)} | while read sample; do
                if grep -v "^#" "bracken/${{{{sample}}}}_kraken2_bracken_report.txt" | grep -q "[0-9]"; then
                    cp "bracken/${{{{sample}}}}_kraken2_bracken_report.txt" "bracken/${{{{sample}}}}_kraken2_bracken_report.tmp"
                    
                    # Fix for kreport2krona.py compatibility
                    echo -e '0\\t0\\t0\\tX\\t0\\t  X' >> "bracken/${{{{sample}}}}_kraken2_bracken_report.tmp"
                    
                    # Convert to Krona format
                    kreport2krona.py \
                        -r "bracken/${{{{sample}}}}_kraken2_bracken_report.tmp" \
                        -o "{{params.outdir}}/${{{{sample}}}}.txt" \
                        2>> {{log}}
                    
                    # Clean taxonomy prefixes
                    sed -i.tmp -E 's/[a-z]__//g' "{{params.outdir}}/${{{{sample}}}}.txt" 2>> {{log}}
                    
                    rm "bracken/${{{{sample}}}}_kraken2_bracken_report.tmp"
                else
                    echo '[{MODULE}] ⚠ Skipping '${{{{sample}}}}' - no Bracken results' | tee -a {{log}}
                fi
            done
            
            # Clean up any .tmp files
            rm -f {{params.outdir}}/*.tmp
            
            # Check if we have any valid Krona input files
            krona_files=$(ls {{params.outdir}}/*.txt 2>/dev/null | wc -l)
            
            if [ "$krona_files" -eq 0 ]; then
                echo '[{MODULE}] ⚠ No valid Krona input files - creating placeholder' | tee -a {{log}}
                echo "<html><body><h1>No valid taxonomic classification results</h1></body></html>" > {{output.krona_html}}
            else
                # Generate Krona plot
                ktImportText {{params.outdir}}/*.txt \
                    -o {{output.krona_html}} \
                    2>> {{log}}
                
                echo '[{MODULE}] ✓ Krona plot generation complete' | tee -a {{log}}
            fi
            '''



# ================================================================================
# MULTIQC REPORT
# ================================================================================


rule multiqc_report:
    """
    Generates comprehensive MultiQC report combining all QC metrics.
    Aggregates results from FastQC, fastp, HISAT2, featureCounts, Kraken2, etc.
    """
    input:
        RUN_RULES
    output:
        html_report = "multiqc/multiqc_report.html"
    params:
        dirs = ['.'],
        outdir="multiqc"
    log: "logs/" + TIMESTAMP + "_multiqc_report.log"
    threads: 1
    shell:
        f'''
        echo '[{MODULE}] Generating MultiQC report...' > {{log}}

        multiqc \
            --force \
            -o {{params.outdir}} \
            {{params.dirs}} \
            &>> {{log}}

        echo '[{MODULE}] ✓ MultiQC report generation complete' >> {{log}}
        '''



# ================================================================================
# WORKFLOW COMPLETION
# ================================================================================


onsuccess:
    # Build list of completed analyses
    analyses = ["FastQC (raw)", "Read preprocessing (fastp)", "FastQC (preprocessed)"]
    if RUNQC_BIOTYPE_COMPOSITION:
        analyses.append("Biotype composition (RSeQC/featureCounts)")
    if RUNQC_SPECIES_COMPOSITION:
        analyses.append("Species composition (Kraken2/Bracken)")
    
    analyses_str = "\n".join([f"[{MODULE}]   • {a}" for a in analyses])
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✓ Quality control workflow completed!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Project ID: {PROJID}'
        echo '[{MODULE}] Processed samples: {len(NAMES)}'
        echo '[{MODULE}] Sequencing type: {'Paired-end' if PAIRED else 'Single-end'}'
        echo '[{MODULE}] Sampling: {'Full dataset' if config['runqc']['number'] == 1 else f"{config['runqc']['number']}{'%' if config['runqc']['number'] < 1 else ' reads'}"}'
        echo ''
        echo '[{MODULE}] Completed analyses:'
        echo '{analyses_str}'
        echo ''
        echo '[{MODULE}] Output directories:'
        echo '[{MODULE}]   - fq_raw_sampled/ (subsampled raw reads)'
        echo '[{MODULE}]   - fastqc_raw/ (raw read QC reports)'
        echo '[{MODULE}]   - fq_preprocessed/ (quality-filtered reads)'
        echo '[{MODULE}]   - fastqc_preprocessed/ (filtered read QC reports)'
        if [ -n "{RUNQC_BIOTYPE_COMPOSITION}" ]; then
            echo '[{MODULE}]   - rseqc/ (strandedness and read distribution)'
            echo '[{MODULE}]   - hisat2_aligned/ (alignment files)'
            echo '[{MODULE}]   - hisat2_aligned_featurecounts/ (biotype counts and plots)'
        fi
        if [ -n "{RUNQC_SPECIES_COMPOSITION}" ]; then
            echo '[{MODULE}]   - kraken2/ (taxonomic classification)'
            echo '[{MODULE}]   - bracken/ (species abundance)'
            echo '[{MODULE}]   - krona/ (interactive taxonomy visualization)'
        fi
        echo '[{MODULE}]   - multiqc/ (comprehensive QC report)'
        echo ''
        echo '[{MODULE}] View comprehensive results: multiqc/multiqc_report.html'
    """)


onerror:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✗ Quality control workflow FAILED!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Common issues:'
        echo '[{MODULE}]   - Insufficient disk space for intermediate files'
        echo '[{MODULE}]   - Memory errors during alignment or classification'
        echo '[{MODULE}]   - Missing or corrupted FASTQ input files'
        echo '[{MODULE}]   - Incorrect reference paths (HISAT2 index, GTF, Kraken2 DB)'
        echo '[{MODULE}]   - Incompatible GTF format for biotype analysis'
        echo '[{MODULE}]   - Missing Bracken database for selected read length'
        echo '[{MODULE}]   - Required tools not in PATH (seqkit, fastp, samtools, etc.)'
        echo ''
        echo '[{MODULE}] Troubleshooting steps:'
        echo '[{MODULE}]   1. Check logs/ directory for detailed error messages'
        echo '[{MODULE}]   2. Verify all input files exist and are readable'
        echo '[{MODULE}]   3. Confirm reference databases are properly indexed'
        echo '[{MODULE}]   4. Ensure sufficient resources (CPU, RAM, disk space)'
        echo ''
        echo '[{MODULE}] Failed at timestamp: {TIMESTAMP}'
        echo ''
    """)
