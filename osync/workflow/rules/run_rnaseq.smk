"""
================================================================================
RUN_RNASEQ MODULE - RNA-seq Analysis Workflow
================================================================================
Performs comprehensive RNA-seq analysis including:
- Read preprocessing with fastp and optional rRNA/contaminant removal (SortMeRNA)
- Alignment to reference genome (HISAT2)
- Gene-level quantification (featureCounts)
- Optional transcript-level quantification (kallisto)
- Optional exon-level quantification with flattened GTF
- Strand-specific BAM splitting and coverage tracks
- Genome coverage visualization (bedGraph format with RPM normalization)
================================================================================
"""

import os
import re
import glob
import subprocess
from pathlib import Path

# Setup workflow environment
workflow_dir = Path(workflow.snakefile).parent.absolute()
sys.path.insert(0, str(workflow_dir.parent))

# Import utilities
from utils.workflow_utils import get_timestamp

# Module identifier
MODULE = 'run_rnaseq'
TIMESTAMP = get_timestamp()


# ================================================================================
# CONFIGURATION AND SETUP
# ================================================================================

# Set variables (quiet parsing - details shown in onstart)
PROJID = os.path.basename(config['run_rnaseq']['indir'])
FILES = glob.glob(os.path.join(config['run_rnaseq']['indir'], '*.f*q.gz'))
EXT = re.search(r'[-\w]+(.f(ast)?q.gz)', FILES[0]).group(1)

try:
    NAMES = list(set([re.search(r'([-\w]+)(?=_[12])_[12].f(ast)?q.gz', x).group(1) for x in FILES]))
except AttributeError:
    NAMES = list(set([re.search(r'([-\w]+).f(ast)?q.gz', x).group(1) for x in FILES]))

PAIRED = False if len(NAMES) == len(FILES) else True

# Set cores and opts
AVAIL_THREADS = config['run_rnaseq']['jobs']

# Parse fastp options
try:
    OPTS = config['run_rnaseq']['fastp_opt']
    FASTP_THREADS = int(re.search(r'(-w|--thread) (\d+)', OPTS).group(2))
    config['run_rnaseq']['fastp_opt'] = re.sub(r'(-w|--thread) (\d+) ','', OPTS)
except AttributeError:
    FASTP_THREADS = 3
FASTP_THREADS = AVAIL_THREADS if FASTP_THREADS > AVAIL_THREADS else FASTP_THREADS

# Parse HISAT2 options
try:
    OPTS = config['run_rnaseq']['hisat2_opt']
    HISAT2_THREADS = int(re.search(r'(-p|--thread) (\d+)', OPTS).group(2))
    config['run_rnaseq']['hisat2_opt'] = re.sub(r'(-p|--thread) (\d+) ', '', OPTS)
except AttributeError:
    HISAT2_THREADS = 3
HISAT2_THREADS = AVAIL_THREADS if HISAT2_THREADS > AVAIL_THREADS else HISAT2_THREADS

# Detect strandedness from HISAT2 options
try:
    OPTS = config['run_rnaseq']['hisat2_opt']
    STRANDED = re.search(r'--rna-strandness(?= *(F|R)) *([FR][FR]?)', OPTS).group(2)
except AttributeError:
    STRANDED = False

# Parse kallisto options if enabled
if config["run_rnaseq"]["kallisto_ref"]:
    try:
        OPTS = config['run_rnaseq']['kallisto_opt']
        KALLISTO_THREADS = int(re.search(r'(-t|--threads) (\d+)', OPTS).group(2))
        config['run_rnaseq']['kallisto_opt'] = re.sub(r'(-t|--threads) (\d+) ', '', OPTS)
    except AttributeError:
        KALLISTO_THREADS = 3
    KALLISTO_THREADS = AVAIL_THREADS if KALLISTO_THREADS > AVAIL_THREADS else KALLISTO_THREADS

    try:
        OPTS = config['run_rnaseq']['kallisto_opt']
        KALLISTO_FRAGMENT_LENGTH = int(re.search(r'(-l|--fragment-length) (\d+)', OPTS).group(2))
        config['run_rnaseq']['kallisto_opt'] = re.sub(r'(-l|--fragment-length) (\d+) ', '', OPTS)
    except AttributeError:
        KALLISTO_FRAGMENT_LENGTH = 100

    try:
        OPTS = config['run_rnaseq']['kallisto_opt']
        KALLISTO_FRAGMENT_LENGTH_SD = int(re.search(r'(-s|--sd) (\d+)', OPTS).group(2))
        config['run_rnaseq']['kallisto_opt'] = re.sub(r'(-s|--sd) (\d+) ', '', OPTS)
    except AttributeError:
        KALLISTO_FRAGMENT_LENGTH_SD = 30

# Check for exon quantification
if config["run_rnaseq"]["exon_quant"]:
    EXON_QUANT = True
else:
    EXON_QUANT = False


# ================================================================================
# WORKFLOW INITIALIZATION
# ================================================================================

onstart:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] RNA-seq Analysis Workflow'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Project: {PROJID}'
        echo '[{MODULE}] Samples: {len(NAMES)}'
        echo '[{MODULE}] Layout: {'Paired-end' if PAIRED else 'Single-end'}'
        echo '[{MODULE}] Strandedness: {STRANDED if STRANDED else 'unstranded'}'
        echo ''
        echo '[{MODULE}] Enabled analyses:'
        echo '[{MODULE}]   • Read preprocessing (fastp)'
        if [ -n "{config['run_rnaseq'].get('removal_ref', '')}" ]; then
            echo '[{MODULE}]   • rRNA/contaminant removal (SortMeRNA)'
        fi
        echo '[{MODULE}]   • Genome alignment (HISAT2)'
        echo '[{MODULE}]   • Gene quantification (featureCounts)'
        if [ -n "{STRANDED}" ]; then
            echo '[{MODULE}]   • Strand-specific BAM splitting'
            echo '[{MODULE}]   • Strand-specific coverage tracks'
        fi
        if [ -n "{config['run_rnaseq'].get('kallisto_ref', '')}" ]; then
            echo '[{MODULE}]   • Transcript quantification (kallisto)'
        fi
        if [ "{config['run_rnaseq'].get('exon_quant', False)}" = "True" ]; then
            echo '[{MODULE}]   • Exon-level quantification'
        fi
        echo '[{MODULE}]   • Genome coverage tracks (bedGraph)'
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
    ["fq_preprocessed/{name}_fastp.json",
    "hisat2_aligned/{name}.bam",
    "hisat2_aligned_featurecounts/" + PROJID + "_featurecounts.txt",
    "hisat2_aligned_bedgraph/{name}.bdg"],
    name=NAMES
)

if STRANDED:
    RUN_RNASEQ_STRANDED = expand(
        ["hisat2_aligned_split/{name}_fwd.bam",
        "hisat2_aligned_split/{name}_rev.bam",
        "hisat2_aligned_split_bedgraph/{name}_fwd.bdg",
        "hisat2_aligned_split_bedgraph/{name}_rev.bdg"],
        name=NAMES
    )
    RUN_RULES.append(RUN_RNASEQ_STRANDED)

RUN_RNASEQ_KALLISTO = None

if config["run_rnaseq"]["kallisto_ref"]:
    RUN_RNASEQ_KALLISTO = expand(
        ["kallisto_transcript_abundance/{name}_abundance.h5"],
        name=NAMES
    )
    RUN_RULES.append(RUN_RNASEQ_KALLISTO)

if config["run_rnaseq"]["exon_quant"]:
    RUN_RNASEQ_EXON_QUANT = [
        "hisat2_aligned_featurecounts_exon_quant/flattened_gtf_merge.saf",
        "hisat2_aligned_featurecounts_exon_quant/flattened_gtf_chop.saf",
        "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_merge.featureCounts",
        "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_chop.featureCounts",
        "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_merge.featureCounts.jcounts",
        "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_chop.featureCounts.jcounts"
    ]
    RUN_RULES.append(RUN_RNASEQ_EXON_QUANT)


rule all:
    input:
        RUN_RULES


# ================================================================================
# READ PREPROCESSING
# ================================================================================

if PAIRED:

    rule preprocess_fq:
        """
        Performs preprocessing on paired-end raw fastq files using fastp.
        Includes adapter trimming, quality filtering, and length filtering.
        Optionally removes rRNA/contaminant reads using SortMeRNA.
        """
        input:
            fq1 = config['run_rnaseq']['indir'] + "/{name}_1" + EXT,
            fq2 = config['run_rnaseq']['indir'] + "/{name}_2" + EXT
        output:
            fq1_preproc = "fq_preprocessed/{name}_preproc_1.fastq.gz",
            fq2_preproc = "fq_preprocessed/{name}_preproc_2.fastq.gz",
            json = "fq_preprocessed/{name}_fastp.json",
            html = "fq_preprocessed/{name}_fastp.html"
        params:
            fastp_opt = config['run_rnaseq']['fastp_opt'],
            removal_ref = config['run_rnaseq']['removal_ref'],
            outdir = "fq_preprocessed"
        log: "logs/" + TIMESTAMP + "_{name}_preprocess_fq.log"
        threads: FASTP_THREADS
        run:
            shell(
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
            )
            if params['removal_ref']:
                shell(
                    f'''
                    echo '[{MODULE}] Removing rRNA/contaminant reads from {{wildcards.name}}...' >> {{log}}

                    # Clean up any previous runs
                    rm -rf {{params.outdir}}/sortmerna/{{wildcards.name}}*
                    rm -f {{params.outdir}}/sortmerna/out/{{wildcards.name}}*

                    sortmerna \
                        --threads {{threads}} \
                        --ref {{params.removal_ref}} \
                        --reads {{output.fq1_preproc}} \
                        --reads {{output.fq2_preproc}} \
                        --workdir {{params.outdir}}/sortmerna \
                        --kvdb {{params.outdir}}/sortmerna/{{wildcards.name}}_kvdb \
                        --readb {{params.outdir}}/sortmerna/{{wildcards.name}}_readb \
                        --fastx \
                        --aligned {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed \
                        --other {{params.outdir}}/{{wildcards.name}} \
                        --paired_in \
                        --out2 \
                        &>> {{log}}

                    # Rename output files to standard format
                    mv {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed_fwd.fq.gz \
                        {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed_1.fastq.gz
                    mv {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed_rev.fq.gz \
                        {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed_2.fastq.gz

                    mv {{params.outdir}}/{{wildcards.name}}_fwd.fq.gz {{output.fq1_preproc}}
                    mv {{params.outdir}}/{{wildcards.name}}_rev.fq.gz {{output.fq2_preproc}}

                    mv {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed.log \
                        {{params.outdir}}/sortmerna/out/{{wildcards.name}}_sortmerna.log

                    echo '[{MODULE}] ✓ rRNA/contaminant removal complete for {{wildcards.name}}' >> {{log}}
                    '''
                )

else:

    rule preprocess_fq:
        """
        Performs preprocessing on single-end raw fastq files using fastp.
        Includes adapter trimming, quality filtering, and length filtering.
        Optionally removes rRNA/contaminant reads using SortMeRNA.
        """
        input:
            fq = config['run_rnaseq']['indir'] + "/{name}" + EXT
        output:
            fq_preproc = "fq_preprocessed/{name}_preproc.fastq.gz",
            json = "fq_preprocessed/{name}_fastp.json",
            html = "fq_preprocessed/{name}_fastp.html"
        params:
            fastp_opt = config['run_rnaseq']['fastp_opt'],
            removal_ref = config['run_rnaseq']['removal_ref'],
            outdir = "fq_preprocessed"
        log: "logs/" + TIMESTAMP + "_{name}_preprocess_fq.log"
        threads: FASTP_THREADS
        run:
            shell(
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
            )
            if params['removal_ref']:
                shell(
                    f'''
                    echo '[{MODULE}] Removing rRNA/contaminant reads from {{wildcards.name}}...' >> {{log}}
                    
                    # Clean up any previous runs
                    rm -rf {{params.outdir}}/sortmerna/{{wildcards.name}}*
                    rm -f {{params.outdir}}/sortmerna/out/{{wildcards.name}}*

                    sortmerna \
                        --threads {{threads}} \
                        --ref {{params.removal_ref}} \
                        --reads {{output.fq_preproc}} \
                        --workdir {{params.outdir}}/sortmerna \
                        --kvdb {{params.outdir}}/sortmerna/{{wildcards.name}}_kvdb \
                        --readb {{params.outdir}}/sortmerna/{{wildcards.name}}_readb \
                        --fastx \
                        --aligned {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed \
                        --other {{params.outdir}}/{{wildcards.name}} \
                        &>> {{log}}

                    # Rename output files to standard format
                    mv {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed.fq.gz \
                        {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed.fastq.gz

                    mv {{params.outdir}}/{{wildcards.name}}.fq.gz {{output.fq_preproc}}

                    mv {{params.outdir}}/sortmerna/out/{{wildcards.name}}_removed.log \
                        {{params.outdir}}/sortmerna/out/{{wildcards.name}}_sortmerna.log

                    echo '[{MODULE}] ✓ rRNA/contaminant removal complete for {{wildcards.name}}' >> {{log}}
                    '''
                )


# ================================================================================
# ALIGNMENT WITH HISAT2
# ================================================================================

if PAIRED:

    rule hisat2_align:
        """
        Aligns paired-end reads to reference genome using HISAT2.
        Generates sorted BAM files with index and alignment statistics.
        """
        input:
            fq1 = "fq_preprocessed/{name}_preproc_1.fastq.gz",
            fq2 = "fq_preprocessed/{name}_preproc_2.fastq.gz"
        output:
            sam = temp("hisat2_aligned/{name}.sam"),
            bam = "hisat2_aligned/{name}.bam",
            bam_idx = "hisat2_aligned/{name}.bam.bai",
            summary = "hisat2_aligned/{name}_hisat2_summary.txt",
            stats = "hisat2_aligned/{name}_stats.txt"
        params:
            hisat2_ref = config['run_rnaseq']['hisat2_ref'],
            hisat2_opt = config['run_rnaseq']['hisat2_opt']
        log: "logs/" + TIMESTAMP + "_{name}_hisat2_align.log"
        threads: HISAT2_THREADS
        shell:
            f'''
            echo '[{MODULE}] Aligning {{wildcards.name}} with HISAT2...' > {{log}}

            hisat2 \
                --threads {{threads}} \
                {{params.hisat2_opt}} \
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

    if STRANDED:

        rule hisat2_split_bam:
            """
            Splits paired-end BAM into forward and reverse strand-specific files.
            Uses SAM flags to separate reads by strand orientation.
            FR: first-read forward, second-read reverse
            RF: first-read reverse, second-read forward
            """
            input:
                bam = "hisat2_aligned/{name}.bam"
            output:
                bam_fwd = "hisat2_aligned_split/{name}_fwd.bam",
                bam_fwd_idx = "hisat2_aligned_split/{name}_fwd.bam.bai",
                bam_rev = "hisat2_aligned_split/{name}_rev.bam",
                bam_rev_idx = "hisat2_aligned_split/{name}_rev.bam.bai"
            params:
                outdir = "hisat2_aligned_split"
            log: "logs/" + TIMESTAMP + "_{name}_hisat2_split_bam.log"
            threads: HISAT2_THREADS
            run:
                if STRANDED == 'FR':
                    pair1 = "fwd"
                    pair2 = "rev"
                elif STRANDED == 'RF':
                    pair1 = "rev"
                    pair2 = "fwd"

                shell(
                    f'''
                    echo '[{MODULE}] Splitting BAM by strand for {{wildcards.name}}...' > {{log}}

                    # First strand: extract reads based on flag combinations
                    # Include first in pair (0x40) and remove read reverse strand (0x10)
                    samtools view \
                        -@ {{threads}} \
                        -b \
                        -f 64 \
                        -F 16 \
                        {{input.bam}} \
                        > {{params.outdir}}/{{wildcards.name}}_{pair1}1.bam.tmp \
                        2>> {{log}}

                    # Include second in pair (0x80) and read reverse strand (0x10)
                    samtools view \
                        -@ {{threads}} \
                        -b \
                        -f 144 \
                        {{input.bam}} \
                        > {{params.outdir}}/{{wildcards.name}}_{pair1}2.bam.tmp \
                        2>> {{log}}

                    samtools merge \
                        -@ {{threads}} \
                        -f \
                        {{params.outdir}}/{{wildcards.name}}_{pair1}1.bam.tmp \
                        {{params.outdir}}/{{wildcards.name}}_{pair1}2.bam.tmp \
                        -o {{params.outdir}}/{{wildcards.name}}_{pair1}.bam.tmp \
                        &>> {{log}}

                    samtools sort \
                        -@ {{threads}} \
                        {{params.outdir}}/{{wildcards.name}}_{pair1}.bam.tmp \
                        -o {{params.outdir}}/{{wildcards.name}}_{pair1}.bam \
                        &>> {{log}}

                    samtools index \
                        {{params.outdir}}/{{wildcards.name}}_{pair1}.bam \
                        &>> {{log}}

                    # Second strand: extract complementary reads
                    # Include second in pair (0x80) and remove read reverse strand (0x10)
                    samtools view \
                        -@ {{threads}} \
                        -b \
                        -f 128 \
                        -F 16 \
                        {{input.bam}} \
                        > {{params.outdir}}/{{wildcards.name}}_{pair2}1.bam.tmp \
                        2>> {{log}}

                    # Include first in pair (0x40) and read reverse strand (0x10) 
                    samtools view \
                        -@ {{threads}} \
                        -b \
                        -f 80 \
                        {{input.bam}} \
                        > {{params.outdir}}/{{wildcards.name}}_{pair2}2.bam.tmp \
                        2>> {{log}}

                    samtools merge \
                        -@ {{threads}} \
                        -f \
                        {{params.outdir}}/{{wildcards.name}}_{pair2}1.bam.tmp \
                        {{params.outdir}}/{{wildcards.name}}_{pair2}2.bam.tmp \
                        -o {{params.outdir}}/{{wildcards.name}}_{pair2}.bam.tmp \
                        &>> {{log}}

                    samtools sort \
                        -@ {{threads}} \
                        {{params.outdir}}/{{wildcards.name}}_{pair2}.bam.tmp \
                        -o {{params.outdir}}/{{wildcards.name}}_{pair2}.bam \
                        &>> {{log}}
                    
                    samtools index \
                        {{params.outdir}}/{{wildcards.name}}_{pair2}.bam \
                        &>> {{log}}

                    # Cleanup temporary files
                    rm {{params.outdir}}/{{wildcards.name}}*.tmp

                    echo '[{MODULE}] ✓ Strand splitting complete for {{wildcards.name}}' >> {{log}}
                    '''
                )
    
    if RUN_RNASEQ_KALLISTO:

        rule kallisto_quant:
            """
            Quantifies transcript abundances from paired-end reads using kallisto.
            Generates H5 and TSV output files with TPM and estimated counts.
            """
            input:
                fq1 = "fq_preprocessed/{name}_preproc_1.fastq.gz",
                fq2 = "fq_preprocessed/{name}_preproc_2.fastq.gz",
                kallisto_ref = config['run_rnaseq']['kallisto_ref']
            output:
                h5 = "kallisto_transcript_abundance/{name}_abundance.h5",
                tsv = "kallisto_transcript_abundance/{name}_abundance.tsv",
                json = "kallisto_transcript_abundance/{name}_run_info.json",
            params:
                kallisto_opt = config['run_rnaseq']['kallisto_opt'],
                outdir = "kallisto_transcript_abundance"
            log: "logs/" + TIMESTAMP + "_{name}_kallisto_quant.log"
            threads: KALLISTO_THREADS
            shell:
                f'''
                echo '[{MODULE}] Quantifying transcripts for {{wildcards.name}} with kallisto...' > {{log}}

                kallisto quant \
                    --threads {{threads}} \
                    {{params.kallisto_opt}} \
                    --index {{input.kallisto_ref}} \
                    --output-dir {{params.outdir}}/{{wildcards.name}} \
                    {{input.fq1}} \
                    {{input.fq2}} \
                    &>> {{log}}
                
                # Move output files to final location
                mv {{params.outdir}}/{{wildcards.name}}/abundance.h5 {{output.h5}}
                mv {{params.outdir}}/{{wildcards.name}}/abundance.tsv {{output.tsv}}
                mv {{params.outdir}}/{{wildcards.name}}/run_info.json {{output.json}}

                rm -r {{params.outdir}}/{{wildcards.name}}

                echo '[{MODULE}] ✓ Kallisto quantification complete for {{wildcards.name}}' >> {{log}}
                '''

else:

    rule hisat2_align:
        """
        Aligns single-end reads to reference genome using HISAT2.
        Generates sorted BAM files with index and alignment statistics.
        """
        input:
            fq = "fq_preprocessed/{name}_preproc.fastq.gz"
        output:
            sam = temp("hisat2_aligned/{name}.sam"),
            bam = "hisat2_aligned/{name}.bam",
            bam_idx = "hisat2_aligned/{name}.bam.bai",
            summary = "hisat2_aligned/{name}_hisat2_summary.txt",
            stats = "hisat2_aligned/{name}_stats.txt"
        params:
            hisat2_ref = config['run_rnaseq']['hisat2_ref'],
            hisat2_opt = config['run_rnaseq']['hisat2_opt']
        log: "logs/" + TIMESTAMP + "_{name}_hisat2_align.log"
        threads: HISAT2_THREADS
        shell:
            f'''
            echo '[{MODULE}] Aligning {{wildcards.name}} with HISAT2...' > {{log}}

            hisat2 \
                --threads {{threads}} \
                {{params.hisat2_opt}} \
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
    
    if STRANDED:

        rule hisat2_split_bam:
            """
            Splits single-end BAM into forward and reverse strand-specific files.
            Uses SAM flag 0x10 (read reverse strand) to separate by orientation.
            """
            input:
                bam = "hisat2_aligned/{name}.bam"
            output:
                bam_fwd = "hisat2_aligned_split/{name}_fwd.bam",
                bam_fwd_idx = "hisat2_aligned_split/{name}_fwd.bam.bai",
                bam_rev = "hisat2_aligned_split/{name}_rev.bam",
                bam_rev_idx = "hisat2_aligned_split/{name}_rev.bam.bai"
            params:
                outdir = "hisat2_aligned_split"
            log: "logs/" + TIMESTAMP + "_{name}_hisat2_split_bam.log"
            threads: HISAT2_THREADS
            run:
                if STRANDED == 'F':
                    read1 = "fwd"
                    read2 = "rev"
                elif STRANDED == 'R':
                    read1 = "rev"
                    read2 = "fwd"

                shell(
                    f'''
                    echo '[{MODULE}] Splitting BAM by strand for {{wildcards.name}}...' > {{log}}

                    # Forward strand: remove read reverse strand (0x10)
                    samtools view \
                        -@ {{threads}} \
                        -b \
                        -F 16 \
                        {{input.bam}} \
                        > {{params.outdir}}/{{wildcards.name}}_{read1}.bam \
                        2>> {{log}}

                    samtools index \
                        {{params.outdir}}/{{wildcards.name}}_{read1}.bam \
                        &>> {{log}}

                    # Reverse strand: include read reverse strand (0x10)
                    samtools view \
                        -@ {{threads}} \
                        -b \
                        -f 16 \
                        {{input.bam}} \
                        > {{params.outdir}}/{{wildcards.name}}_{read2}.bam \
                        2>> {{log}}
                    
                    samtools index \
                        {{params.outdir}}/{{wildcards.name}}_{read2}.bam \
                        &>> {{log}}

                    echo '[{MODULE}] ✓ Strand splitting complete for {{wildcards.name}}' >> {{log}}
                    '''
                )

    if RUN_RNASEQ_KALLISTO:

        rule kallisto_quant:
            """
            Quantifies transcript abundances from single-end reads using kallisto.
            Generates H5 and TSV output files with TPM and estimated counts.
            Requires fragment length and standard deviation parameters for SE mode.
            """
            input:
                fq = "fq_preprocessed/{name}_preproc.fastq.gz",
                kallisto_ref = config['run_rnaseq']['kallisto_ref']
            output:
                h5 = "kallisto_transcript_abundance/{name}_abundance.h5",
                tsv = "kallisto_transcript_abundance/{name}_abundance.tsv",
                json = "kallisto_transcript_abundance/{name}_run_info.json"
            params:
                kallisto_opt = config['run_rnaseq']['kallisto_opt'],
                frag_len = KALLISTO_FRAGMENT_LENGTH,
                frag_len_sd = KALLISTO_FRAGMENT_LENGTH_SD,
                outdir = "kallisto_transcript_abundance"
            log: "logs/" + TIMESTAMP + "_{name}_kallisto_quant.log"
            threads: KALLISTO_THREADS
            shell:
                f'''
                echo '[{MODULE}] Quantifying transcripts for {{wildcards.name}} with kallisto...' > {{log}}

                kallisto quant \
                    --threads {{threads}} \
                    {{params.kallisto_opt}} \
                    --single \
                    --fragment-length {{params.frag_len}} \
                    --sd {{params.frag_len_sd}} \
                    --index {{input.kallisto_ref}} \
                    --output-dir {{params.outdir}}/{{wildcards.name}} \
                    {{input.fq}} \
                    &>> {{log}}
                
                # Move output files to final location
                mv {{params.outdir}}/{{wildcards.name}}/abundance.h5 {{output.h5}}
                mv {{params.outdir}}/{{wildcards.name}}/abundance.tsv {{output.tsv}}
                mv {{params.outdir}}/{{wildcards.name}}/run_info.json {{output.json}}

                rm -r {{params.outdir}}/{{wildcards.name}}

                echo '[{MODULE}] ✓ Kallisto quantification complete for {{wildcards.name}}' >> {{log}}
                '''


# ================================================================================
# GENE-LEVEL QUANTIFICATION
# ================================================================================

rule feature_counts:
    """
    Gene expression quantification with featureCounts (Subread package).
    Counts reads overlapping gene exons using GTF annotation.
    Automatically detects strandedness and paired-end mode.
    """
    input:
        bams = expand(["hisat2_aligned/{name}.bam"], name=NAMES),
        gtf = config['run_rnaseq']['gtf']
    output:
        featurecounts = "hisat2_aligned_featurecounts/" + PROJID + "_featurecounts.txt"
    params:
        outdir = "hisat2_aligned_featurecounts"
    log: "logs/" + TIMESTAMP + "_feature_counts.log"
    threads: AVAIL_THREADS
    run:
        pe = '-p --countReadPairs -C -B' if PAIRED else ''

        featurecounts_strandedness = '-s 0'
        if STRANDED in ['F','FR']:
            featurecounts_strandedness = '-s 1'
        elif STRANDED in ['R','RF']:
            featurecounts_strandedness = '-s 2'

        shell(
            f'''
            echo '[{MODULE}] Quantifying gene expression with featureCounts...' > {{log}}

            featureCounts \
                -T {{threads}} \
                {pe} \
                {featurecounts_strandedness} \
                -t exon \
                -g gene_id \
                -a {{input.gtf}} \
                -o {{output.featurecounts}} \
                {{input.bams}} \
                &>> {{log}}

            echo '[{MODULE}] ✓ Gene expression quantification complete' >> {{log}}
            '''
        )


# ================================================================================
# EXON-LEVEL QUANTIFICATION (OPTIONAL)
# ================================================================================

if EXON_QUANT:
    rule flatten_gtf:
        """
        Flattens GTF annotation for robust exon-level quantification.
        Generates two SAF formats:
        - merge: Overlapping exons merged into single features
        - chop: Overlapping exons split at boundaries (disjoint bins)
        """
        input:
            gtf = config['run_rnaseq']['gtf']
        output:
            saf_merge = "hisat2_aligned_featurecounts_exon_quant/flattened_gtf_merge.saf",
            saf_chop = "hisat2_aligned_featurecounts_exon_quant/flattened_gtf_chop.saf"
        params:
            outdir = "hisat2_aligned_featurecounts_exon_quant"
        log: "logs/" + TIMESTAMP + "_flatten_gtf.log"
        threads: 1
        shell:
            f'''
            echo '[{MODULE}] Flattening GTF annotation...' > {{log}}

            flattenGTF -a {{input.gtf}} -o {{output.saf_merge}} -t exon -g gene_id \
                &>> {{log}}
            
            flattenGTF -a {{input.gtf}} -o {{output.saf_chop}} -t exon -g gene_id -C \
                &>> {{log}}

            echo '[{MODULE}] ✓ GTF flattening complete' >> {{log}}
            '''
    
    rule feature_counts_exon:
        """
        Exon-level quantification with featureCounts using flattened annotations.
        Generates four count matrices:
        - merge/chop without junction counts (--nonSplitOnly)
        - merge/chop with junction counts (-J flag)
        Allows multi-overlapping features (-O) for comprehensive quantification.
        """
        input:
            bams = expand(["hisat2_aligned/{name}.bam"], name=NAMES),
            saf_merge = "hisat2_aligned_featurecounts_exon_quant/flattened_gtf_merge.saf",
            saf_chop = "hisat2_aligned_featurecounts_exon_quant/flattened_gtf_chop.saf",
            genome = config['run_rnaseq']['genome']
        output:
            featurecounts_merge = "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_merge.featureCounts",
            featurecounts_chop = "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_chop.featureCounts",
            featurecounts_merge_jcounts = "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_merge.featureCounts.jcounts",
            featurecounts_chop_jcounts = "hisat2_aligned_featurecounts_exon_quant/" + PROJID + "_featurecounts_exon_chop.featureCounts.jcounts"
        params:
            outdir = "hisat2_aligned_featurecounts_exon_quant"
        log: "logs/" + TIMESTAMP + "_feature_counts_exon.log"
        threads: AVAIL_THREADS
        run:
            pe = '-p --countReadPairs -C -B' if PAIRED else ''
            featurecounts_strandedness = '-s 0'
            if STRANDED in ['F','FR']:
                featurecounts_strandedness = '-s 1'
            elif STRANDED in ['R','RF']:
                featurecounts_strandedness = '-s 2'
            g = '-G ' + input.genome if input.genome else ''
            
            shell(
                f'''
                echo '[{MODULE}] Quantifying exon expression with flattened GTF...' > {{log}}
                
                # Merge mode - exon counts only
                featureCounts \
                    -T {{threads}} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {{input.saf_merge}} \
                    -o {{output.featurecounts_merge}} \
                    --nonSplitOnly \
                    {{input.bams}} \
                    &>> {{log}}
                
                # Chop mode - exon counts only
                featureCounts \
                    -T {{threads}} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {{input.saf_chop}} \
                    -o {{output.featurecounts_chop}} \
                    --nonSplitOnly \
                    {{input.bams}} \
                    &>> {{log}}
                
                # Merge mode - with junction counts
                featureCounts \
                    -T {{threads}} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {{input.saf_merge}} \
                    -o {{output.featurecounts_merge}} \
                    -J \
                    {g} \
                    {{input.bams}} \
                    &>> {{log}}
                
                # Chop mode - with junction counts
                featureCounts \
                    -T {{threads}} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {{input.saf_chop}} \
                    -o {{output.featurecounts_chop}} \
                    -J \
                    {g} \
                    {{input.bams}} \
                    &>> {{log}}
                
                echo '[{MODULE}] ✓ Exon-level quantification complete' >> {{log}}
                '''
            )


# ================================================================================
# GENOME COVERAGE TRACKS
# ================================================================================

rule make_bedgraph:
    """
    Generates RPM-normalized genome coverage tracks in bedGraph format.
    Includes track line for direct loading into genome browsers.
    Splits reads at splice junctions (-split) for accurate coverage.
    """
    input:
        bam = "hisat2_aligned/{name}.bam"
    output:
        bdg = "hisat2_aligned_bedgraph/{name}.bdg"
    log: "logs/" + TIMESTAMP + "_{name}_make_bedgraph.log"
    threads: 1
    run:
        # Calculate RPM scaling factor
        reads = subprocess.check_output(
            f"samtools view -c {input.bam}", 
            shell=True, 
            text=True
        ).strip()
        scale_factor = 1e6/(int(reads)/2) if PAIRED else 1e6/int(reads)

        pe = '-pc' if PAIRED else ''

        shell(
            f'''
            echo '[{MODULE}] Generating bedGraph for {{wildcards.name}}...' > {{log}}

            bedtools genomecov \
                {pe} \
                -bga \
                -split \
                -trackline \
                -trackopts 'name={{wildcards.name}} color=128,128,128' \
                -scale {scale_factor} \
                -ibam {{input.bam}} \
                > {{output.bdg}} \
                2>> {{log}}

            echo '[{MODULE}] ✓ bedGraph generation complete for {{wildcards.name}}' >> {{log}}
            '''
        )


if STRANDED:

    rule make_stranded_bedgraph:
        """
        Generates strand-specific RPM-normalized genome coverage tracks.
        Creates separate bedGraph files for forward (red) and reverse (blue) strands.
        Normalization factor calculated from total mapped reads across both strands.
        """
        input:
            bam_fwd = "hisat2_aligned_split/{name}_fwd.bam",
            bam_rev = "hisat2_aligned_split/{name}_rev.bam"
        output:
            bdg_fwd = "hisat2_aligned_split_bedgraph/{name}_fwd.bdg",
            bdg_rev = "hisat2_aligned_split_bedgraph/{name}_rev.bdg"
        log: "logs/" + TIMESTAMP + "_{name}_make_stranded_bedgraph.log"
        threads: 1
        run:
            # Calculate RPM scaling factor from total reads
            reads_fwd = subprocess.check_output(
                f"samtools view -c {input.bam_fwd}", 
                shell=True, 
                text=True
            ).strip()
            reads_rev = subprocess.check_output(
                f"samtools view -c {input.bam_rev}", 
                shell=True, 
                text=True
            ).strip()
            
            total_reads = int(reads_fwd) + int(reads_rev)
            scale_factor = 1e6/(total_reads/2) if PAIRED else 1e6/total_reads

            pe = '-pc' if PAIRED else ''

            shell(
                f'''
                echo '[{MODULE}] Generating strand-specific bedGraphs for {{wildcards.name}}...' > {{log}}

                # Forward strand
                bedtools genomecov \
                    {pe} \
                    -bga \
                    -split \
                    -trackline \
                    -trackopts 'name={{wildcards.name}}_fwd color=255,0,0' \
                    -scale {scale_factor} \
                    -ibam {{input.bam_fwd}} \
                    > {{output.bdg_fwd}} \
                    2>> {{log}}

                # Reverse strand
                bedtools genomecov \
                    {pe} \
                    -bga \
                    -split \
                    -trackline \
                    -trackopts 'name={{wildcards.name}}_rev color=0,0,255' \
                    -scale {scale_factor} \
                    -ibam {{input.bam_rev}} \
                    > {{output.bdg_rev}} \
                    2>> {{log}}

                echo '[{MODULE}] ✓ Strand-specific bedGraph generation complete for {{wildcards.name}}' >> {{log}}
                '''
            )


# ================================================================================
# WORKFLOW COMPLETION
# ================================================================================

onsuccess:
    # Build list of completed analyses
    analyses = ["Read preprocessing (fastp)", "Genome alignment (HISAT2)", "Gene quantification (featureCounts)", "Genome coverage tracks (bedGraph)"]
    if config['run_rnaseq'].get('removal_ref'):
        analyses.insert(1, "rRNA/contaminant removal (SortMeRNA)")
    if STRANDED:
        analyses.append("Strand-specific BAM splitting")
        analyses.append("Strand-specific coverage tracks")
    if config['run_rnaseq'].get('kallisto_ref'):
        analyses.append("Transcript quantification (kallisto)")
    if config['run_rnaseq'].get('exon_quant'):
        analyses.append("Exon-level quantification")
    
    analyses_str = "\n".join([f"[{MODULE}]   • {a}" for a in analyses])
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✓ RNA-seq analysis workflow completed!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Project ID: {PROJID}'
        echo '[{MODULE}] Processed samples: {len(NAMES)}'
        echo '[{MODULE}] Sequencing type: {'Paired-end' if PAIRED else 'Single-end'}'
        echo '[{MODULE}] Strandedness: {STRANDED if STRANDED else 'unstranded'}'
        echo ''
        echo '[{MODULE}] Completed analyses:'
        echo '{analyses_str}'
        echo ''
        echo '[{MODULE}] Output directories:'
        echo '[{MODULE}]   - fq_preprocessed/ (quality-filtered reads)'
        if [ -n "{config['run_rnaseq'].get('removal_ref', '')}" ]; then
            echo '[{MODULE}]   - fq_preprocessed/sortmerna/out/ (removed rRNA/contaminants)'
        fi
        echo '[{MODULE}]   - hisat2_aligned/ (alignment BAM files)'
        echo '[{MODULE}]   - hisat2_aligned_featurecounts/ (gene expression counts)'
        echo '[{MODULE}]   - hisat2_aligned_bedgraph/ (genome coverage tracks)'
        if [ -n "{STRANDED}" ]; then
            echo '[{MODULE}]   - hisat2_aligned_split/ (strand-specific BAM files)'
            echo '[{MODULE}]   - hisat2_aligned_split_bedgraph/ (strand-specific coverage)'
        fi
        if [ -n "{config['run_rnaseq'].get('kallisto_ref', '')}" ]; then
            echo '[{MODULE}]   - kallisto_transcript_abundance/ (transcript TPMs)'
        fi
        if [ "{config['run_rnaseq'].get('exon_quant', False)}" = "True" ]; then
            echo '[{MODULE}]   - hisat2_aligned_featurecounts_exon_quant/ (exon counts)'
        fi
        echo ''
        echo '[{MODULE}] Primary outputs:'
        echo '[{MODULE}]   - Gene counts: hisat2_aligned_featurecounts/{PROJID}_featurecounts.txt'
        if [ -n "{config['run_rnaseq'].get('kallisto_ref', '')}" ]; then
            echo '[{MODULE}]   - Transcript abundance: kallisto_transcript_abundance/*.tsv'
        fi
        echo '[{MODULE}]   - Coverage tracks: hisat2_aligned_bedgraph/*.bdg'
    """)

onerror:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✗ RNA-seq analysis workflow FAILED!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Common issues:'
        echo '[{MODULE}]   - Insufficient disk space for BAM files'
        echo '[{MODULE}]   - Memory errors during alignment'
        echo '[{MODULE}]   - Missing or corrupted FASTQ input files'
        echo '[{MODULE}]   - Incorrect HISAT2 index or GTF file paths'
        echo '[{MODULE}]   - Incompatible strandedness settings'
        echo '[{MODULE}]   - SortMeRNA database issues'
        echo '[{MODULE}]   - Required tools not in PATH (fastp, hisat2, samtools, etc.)'
        echo ''
        echo '[{MODULE}] Troubleshooting steps:'
        echo '[{MODULE}]   1. Check logs/ directory for detailed error messages'
        echo '[{MODULE}]   2. Verify all input files exist and are readable'
        echo '[{MODULE}]   3. Confirm reference files (index, GTF) are properly formatted'
        echo '[{MODULE}]   4. Ensure sufficient resources (CPU, RAM, disk space)'
        echo '[{MODULE}]   5. Validate strandedness settings match library prep protocol'
        echo ''
        echo '[{MODULE}] Failed at timestamp: {TIMESTAMP}'
        echo ''
    """)
