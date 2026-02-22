"""
================================================================================
RUN_ATACSEQ MODULE - ATAC-seq Analysis Workflow
================================================================================
Performs comprehensive ATAC-seq analysis including:
- Read preprocessing with fastp (adapter trimming, quality filtering)
- Alignment to reference genome (HISAT2, no spliced alignment)
- Post-alignment processing (filtering, optional Tn5 shifting)
- Peak calling with MACS2
- Quantification of reads in peaks (featureCounts)
- Quantification of reads in genes (featureCounts)
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
MODULE = 'run_atacseq'
TIMESTAMP = get_timestamp()


# ================================================================================
# CONFIGURATION AND SETUP
# ================================================================================

# Set variables
PROJID = os.path.basename(config['run_atacseq']['indir'])
FILES = glob.glob(os.path.join(config['run_atacseq']['indir'], '*.f*q.gz'))
EXT = re.search(r'[-\w]+(.f(ast)?q.gz)', FILES[0]).group(1)

try:
    NAMES = list(set([re.search(r'([-\w]+)(?=_[12])_[12].f(ast)?q.gz', x).group(1) for x in FILES]))
except AttributeError:
    NAMES = list(set([re.search(r'([-\w]+).f(ast)?q.gz', x).group(1) for x in FILES]))

PAIRED = False if len(NAMES) == len(FILES) else True

# Set cores and opts
AVAIL_THREADS = config['run_atacseq']['jobs']

# Parse fastp options
try:
    OPTS = config['run_atacseq']['fastp_opt']
    FASTP_THREADS = int(re.search(r'(-w|--thread) (\d+)', OPTS).group(2))
    config['run_atacseq']['fastp_opt'] = re.sub(r'(-w|--thread) (\d+) ','', OPTS)
except AttributeError:
    FASTP_THREADS = 3
FASTP_THREADS = AVAIL_THREADS if FASTP_THREADS > AVAIL_THREADS else FASTP_THREADS

# Parse HISAT2 options
try:
    OPTS = config['run_atacseq']['hisat2_opt']
    HISAT2_THREADS = int(re.search(r'(-p|--thread) (\d+)', OPTS).group(2))
    config['run_atacseq']['hisat2_opt'] = re.sub(r'(-p|--thread) (\d+) ', '', OPTS)
except AttributeError:
    HISAT2_THREADS = 3
HISAT2_THREADS = AVAIL_THREADS if HISAT2_THREADS > AVAIL_THREADS else HISAT2_THREADS


# ================================================================================
# WORKFLOW INITIALIZATION
# ================================================================================

onstart:
    # Build configuration summary
    layout = "Paired-end" if PAIRED else "Single-end"
    tn5_status = "Enabled (+4bp fwd, -5bp rev)" if config['run_atacseq']['tn5_shift'] else "Disabled"
    exclusion_status = "Custom exclusion list + mitochondria" if config['run_atacseq']['exclusion_list'] else "Mitochondria only"
    
    # List sample names
    samples_list = "\n".join([f"[{MODULE}]   • {name}" for name in NAMES])
    
    # Build workflow steps
    workflow_steps = [
        "[{MODULE}]   1. Read preprocessing (fastp)",
        "[{MODULE}]   2. Genome alignment (HISAT2, no splicing)",
        "[{MODULE}]   3. Post-alignment filtering (mitochondria/blacklist)"
    ]
    if config['run_atacseq']['tn5_shift']:
        workflow_steps.append("[{MODULE}]   4. Tn5 transposase shifting")
        workflow_steps.append("[{MODULE}]   5. Peak calling (MACS2)")
        workflow_steps.append("[{MODULE}]   6. Peak quantification (featureCounts)")
        workflow_steps.append("[{MODULE}]   7. Gene quantification (featureCounts)")
        workflow_steps.append("[{MODULE}]   8. Coverage tracks (bedGraph, RPM-normalized)")
    else:
        workflow_steps.append("[{MODULE}]   4. Peak calling (MACS2)")
        workflow_steps.append("[{MODULE}]   5. Peak quantification (featureCounts)")
        workflow_steps.append("[{MODULE}]   6. Gene quantification (featureCounts)")
        workflow_steps.append("[{MODULE}]   7. Coverage tracks (bedGraph, RPM-normalized)")
    
    steps_str = "\n".join(workflow_steps)
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ATAC-seq Analysis Workflow'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Project: {PROJID}'
        echo '[{MODULE}] Samples: {len(NAMES)}'
        echo '[{MODULE}] Layout: {layout}'
        echo ''
        echo '[{MODULE}] Configuration:'
        echo '[{MODULE}]   • Threads (fastp): {FASTP_THREADS}'
        echo '[{MODULE}]   • Threads (HISAT2): {HISAT2_THREADS}'
        echo '[{MODULE}]   • Threads (other): {AVAIL_THREADS}'
        echo '[{MODULE}]   • Tn5 shifting: {tn5_status}'
        echo '[{MODULE}]   • Exclusion regions: {exclusion_status}'
        echo '[{MODULE}]   • Genome size: {config["run_atacseq"]["genome_size"]}'
        echo ''
        echo '[{MODULE}] Sample names:'
        echo '{samples_list}'
        echo ''
        echo '[{MODULE}] Workflow steps:'
        echo '{steps_str}'
        echo ''
        echo '[{MODULE}] Reference files:'
        echo '[{MODULE}]   • HISAT2 index: {os.path.basename(config["run_atacseq"]["hisat2_ref"])}'
        echo '[{MODULE}]   • GTF annotation: {os.path.basename(config["run_atacseq"]["gtf"])}'
        echo '[{MODULE}]   • Chromosome sizes: {os.path.basename(config["run_atacseq"]["chr_sizes"])}'
        if [ -n "{config['run_atacseq'].get('exclusion_list', '')}" ]; then
            echo '[{MODULE}]   • Exclusion list: {os.path.basename(config["run_atacseq"]["exclusion_list"])}'
        fi
        echo ''
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
    "hisat2_aligned_processed/{name}_processed.bam",
    "macs2_peaks/{name}_processed_peaks.narrowPeak",
    "featurecounts_reads_in_peaks/{name}_reads_in_peaks.txt",
    "featurecounts_reads_in_genes/" + PROJID + "_reads_in_genes.txt",
    "scaled_bedgraph/{name}_processed.bdg"],
    name=NAMES
)

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
        ATAC-seq specific: typically removes adapters and low-quality bases.
        """
        input:
            fq1 = config['run_atacseq']['indir'] + "/{name}_1" + EXT,
            fq2 = config['run_atacseq']['indir'] + "/{name}_2" + EXT
        output:
            fq1_preproc = "fq_preprocessed/{name}_preproc_1.fastq.gz",
            fq2_preproc = "fq_preprocessed/{name}_preproc_2.fastq.gz",
            json = "fq_preprocessed/{name}_fastp.json",
            html = "fq_preprocessed/{name}_fastp.html"
        params:
            fastp_opt = config['run_atacseq']['fastp_opt'],
            outdir = "fq_preprocessed"
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
        ATAC-seq specific: typically removes adapters and low-quality bases.
        """
        input:
            fq = config['run_atacseq']['indir'] + "/{name}" + EXT
        output:
            fq_preproc = "fq_preprocessed/{name}_preproc.fastq.gz",
            json = "fq_preprocessed/{name}_fastp.json",
            html = "fq_preprocessed/{name}_fastp.html"
        params:
            fastp_opt = config['run_atacseq']['fastp_opt'],
            outdir = "fq_preprocessed"
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


# ================================================================================
# ALIGNMENT WITH HISAT2
# ================================================================================

if PAIRED:

    rule hisat2_align:
        """
        Aligns paired-end reads to reference genome using HISAT2.
        ATAC-seq specific: disables spliced alignment and temp splicesite detection.
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
            hisat2_ref = config['run_atacseq']['hisat2_ref'],
            hisat2_opt = config['run_atacseq']['hisat2_opt']
        log: "logs/" + TIMESTAMP + "_{name}_hisat2_align.log"
        threads: HISAT2_THREADS
        shell:
            f'''
            echo '[{MODULE}] Aligning {{wildcards.name}} with HISAT2 (no spliced alignment)...' > {{log}}

            hisat2 \
                --threads {{threads}} \
                {{params.hisat2_opt}} \
                -1 {{input.fq1}} \
                -2 {{input.fq2}} \
                -x {{params.hisat2_ref}} \
                --no-spliced-alignment \
                --no-temp-splicesite \
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

else:

    rule hisat2_align:
        """
        Aligns single-end reads to reference genome using HISAT2.
        ATAC-seq specific: disables spliced alignment and temp splicesite detection.
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
            hisat2_ref = config['run_atacseq']['hisat2_ref'],
            hisat2_opt = config['run_atacseq']['hisat2_opt']
        log: "logs/" + TIMESTAMP + "_{name}_hisat2_align.log"
        threads: HISAT2_THREADS
        shell:
            f'''
            echo '[{MODULE}] Aligning {{wildcards.name}} with HISAT2 (no spliced alignment)...' > {{log}}

            hisat2 \
                --threads {{threads}} \
                {{params.hisat2_opt}} \
                -U {{input.fq}} \
                -x {{params.hisat2_ref}} \
                --no-spliced-alignment \
                --no-temp-splicesite \
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


# ================================================================================
# POST-ALIGNMENT PROCESSING
# ================================================================================

rule generate_inclusion_list:
    """
    Generates list of genomic regions to include based on exclusion list.
    ATAC-seq specific: excludes mitochondrial DNA and optionally blacklisted regions.
    Uses complement of exclusion list to create inclusion regions.
    """
    input:
        exclusion_list = config['run_atacseq']['exclusion_list'],
        chr_sizes = config['run_atacseq']['chr_sizes']
    output:
        "hisat2_aligned_processed/inclusion_list.bed",
    params:
        outdir = "hisat2_aligned_processed"
    log: "logs/" + TIMESTAMP + "_generate_inclusion_list.log"
    threads: 1
    run:
        if input.exclusion_list:
            shell(
                f'''
                echo '[{MODULE}] Generating inclusion list from exclusion regions...' > {{log}}

                # Add mitochondrial regions to exclusion
                awk '$1 ~ /(MT|mt|chrMT|chrmt|chrM|chrm)/ {{{{print $1"\\t0\\t"$2}}}}' {{input.chr_sizes}} \
                    > {{params.outdir}}/exclusion_list.tmp \
                    2>> {{log}}
                
                # Add custom exclusion list
                awk '{{{{print $1,$2,$3}}}}' OFS='\\t' {{input.exclusion_list}} \
                    >> {{params.outdir}}/exclusion_list.tmp \
                    2>> {{log}}

                # Sort and generate complement (inclusion list)
                sortBed -g {{input.chr_sizes}} -i {{params.outdir}}/exclusion_list.tmp \
                    > {{params.outdir}}/exclusion_list_sorted.tmp \
                    2>> {{log}}
                
                complementBed -g {{input.chr_sizes}} -i {{params.outdir}}/exclusion_list_sorted.tmp \
                    > {{params.outdir}}/inclusion_list.bed \
                    2>> {{log}}
                
                rm {{params.outdir}}/*.tmp

                echo '[{MODULE}] ✓ Inclusion list generated (excluding mitochondria + custom regions)' >> {{log}}
                '''
            )
        else:
            shell(
                f'''
                echo '[{MODULE}] Generating inclusion list (mitochondria only)...' > {{log}}

                # Add mitochondrial regions to exclusion
                awk '$1 ~ /(MT|mt|chrMT|chrmt|chrM|chrm)/ {{{{print $1"\\t0\\t"$2}}}}' {{input.chr_sizes}} \
                    > {{params.outdir}}/exclusion_list.tmp \
                    2>> {{log}}

                # Sort and generate complement (inclusion list)
                sortBed -g {{input.chr_sizes}} -i {{params.outdir}}/exclusion_list.tmp \
                    > {{params.outdir}}/exclusion_list_sorted.tmp \
                    2>> {{log}}
                
                complementBed -g {{input.chr_sizes}} -i {{params.outdir}}/exclusion_list_sorted.tmp \
                    > {{params.outdir}}/inclusion_list.bed \
                    2>> {{log}}
                
                rm {{params.outdir}}/*.tmp

                echo '[{MODULE}] ✓ Inclusion list generated (excluding mitochondria only)' >> {{log}}
                '''
            )


rule hisat2_align_process:
    """
    Filters aligned reads to inclusion regions and optionally performs Tn5 shifting.
    ATAC-seq specific: 
    - Removes reads in excluded regions (mitochondria, blacklisted)
    - Optional Tn5 transposase shift (+4bp forward, -5bp reverse strand)
    """
    input:
        bam = "hisat2_aligned/{name}.bam",
        inclusion_list = "hisat2_aligned_processed/inclusion_list.bed"
    output:
        processed_bam = "hisat2_aligned_processed/{name}_processed.bam",
    params:
        outdir = "hisat2_aligned_processed"
    log: "logs/" + TIMESTAMP + "_{name}_hisat2_align_process.log"
    threads: HISAT2_THREADS
    run:
        shell(f" echo '[{MODULE}] Processing aligned reads for {{wildcards.name}}...' > {{log}} ")

        # Filter to inclusion regions
        shell(
            f'''
            echo '[{MODULE}] Filtering to inclusion regions...' >> {{log}}

            samtools view \
                -@ {{threads}} \
                -b \
                -L {{input.inclusion_list}} \
                {{input.bam}} \
                > {{output.processed_bam}} \
                2>> {{log}}
            
            samtools sort \
                -@ {{threads}} \
                {{output.processed_bam}} \
                -o {{output.processed_bam}} \
                &>> {{log}}
            
            samtools index \
                -@ {{threads}} \
                {{output.processed_bam}} \
                &>> {{log}}
            '''
        )

        # Optional Tn5 shifting
        if config['run_atacseq']['tn5_shift']:
            shell(
                f'''
                echo '[{MODULE}] Applying Tn5 transposase shift (+4bp fwd, -5bp rev)...' >> {{log}}

                alignmentSieve \
                    -p {{threads}} \
                    --ATACshift \
                    -b {{output.processed_bam}} \
                    -o {{output.processed_bam}} \
                    &>> {{log}}

                samtools sort \
                    -@ {{threads}} \
                    {{output.processed_bam}} \
                    -o {{output.processed_bam}} \
                    &>> {{log}}
                
                samtools index \
                    -@ {{threads}} \
                    {{output.processed_bam}} \
                    &>> {{log}}
                '''
            )
        
        shell(f" echo '[{MODULE}] ✓ Processing complete for {{wildcards.name}}' >> {{log}} ")


# ================================================================================
# PEAK CALLING AND QUANTIFICATION
# ================================================================================

rule macs2_call_peaks:
    """
    Calls accessible chromatin peaks using MACS2.
    ATAC-seq specific: uses appropriate format (BAMPE for paired-end, BAM for single-end).
    Generates narrowPeak files and bedGraph coverage tracks.
    """
    input:
        bam = "hisat2_aligned_processed/{name}_processed.bam"
    output:
        peaks_bed = "macs2_peaks/{name}_processed_peaks.narrowPeak"
    params:
        genome_size = config['run_atacseq']['genome_size'],
        outdir = "macs2_peaks"
    log: "logs/" + TIMESTAMP + "_{name}_macs2_call_peaks.log"
    threads: 1
    run:
        fmt = 'BAMPE' if PAIRED else 'BAM'

        shell(
            f'''
            echo '[{{MODULE}}] Calling peaks with MACS2 (format: {fmt})...' > {{log}}

            macs2 callpeak \
                --bdg \
                -f {fmt} \
                -g {{params.genome_size}} \
                -t {{input.bam}} \
                -n {{wildcards.name}}_processed \
                --outdir {{params.outdir}} \
                &>> {{log}}

            echo '[{MODULE}] ✓ Peak calling complete for {{wildcards.name}}' >> {{log}}
            '''
        )


rule reads_in_peaks:
    """
    Quantifies reads/fragments falling into called peak regions using featureCounts.
    Creates SAF format annotation from narrowPeak file.
    Provides per-sample peak accessibility quantification.
    """
    input:
        bam = "hisat2_aligned_processed/{name}_processed.bam",
        peaks_bed = "macs2_peaks/{name}_processed_peaks.narrowPeak"
    output:
        saf = temp("featurecounts_reads_in_peaks/{name}_peaks.saf"),
        reads_in_peaks = "featurecounts_reads_in_peaks/{name}_reads_in_peaks.txt"
    params:
        outdir = "featurecounts_reads_in_peaks"
    log: "logs/" + TIMESTAMP + "_{name}_reads_in_peaks.log"
    threads: HISAT2_THREADS
    run:
        pe = '-p --countReadPairs -C -B' if PAIRED else ''

        shell(
            f'''
            echo '[{MODULE}] Quantifying reads in peaks for {{wildcards.name}}...' > {{log}}

            # Convert narrowPeak to SAF format
            echo -e 'GeneID\\tChr\\tStart\\tEnd\\tStrand' > {{output.saf}} 2>> {{log}}
            awk '{{{{print $4,$1,$2+1,$3,"."}}}}' OFS='\\t' {{input.peaks_bed}} \
                >> {{output.saf}} \
                2>> {{log}}

            # Count reads in peaks
            featureCounts \
                -T {{threads}} \
                {pe} \
                -F SAF \
                -a {{output.saf}} \
                -o {{output.reads_in_peaks}} \
                {{input.bam}} \
                &>> {{log}}

            echo '[{MODULE}] ✓ Peak quantification complete for {{wildcards.name}}' >> {{log}}
            '''
        )


rule reads_in_genes:
    """
    Quantifies reads/fragments falling into annotated genes using featureCounts.
    Useful for assessing transcription start site (TSS) enrichment.
    Processes all samples together for consistent gene-level quantification.
    """
    input:
        bams = expand(["hisat2_aligned_processed/{name}_processed.bam"], name=NAMES),
        gtf = config['run_atacseq']['gtf']
    output:
        reads_in_genes = "featurecounts_reads_in_genes/" + PROJID + "_reads_in_genes.txt"
    params:
        outdir = "featurecounts_reads_in_genes"
    log: "logs/" + TIMESTAMP + "_reads_in_genes.log"
    threads: AVAIL_THREADS
    run:
        pe = '-p --countReadPairs -C -B' if PAIRED else ''

        shell(
            f'''
            echo '[{MODULE}] Quantifying reads in annotated genes...' > {{log}}

            featureCounts \
                -T {{threads}} \
                {pe} \
                -t gene \
                -g gene_id \
                -a {{input.gtf}} \
                -o {{output.reads_in_genes}} \
                {{input.bams}} \
                &>> {{log}}

            echo '[{MODULE}] ✓ Gene quantification complete for all samples' >> {{log}}
            '''
        )


# ================================================================================
# VISUALIZATION
# ================================================================================

rule make_scaled_bedgraph:
    """
    Generates RPM-normalized genome coverage files in bedGraph format.
    ATAC-seq specific: visualizes open chromatin landscape.
    Scaling factor: 1,000,000 / (total reads / 2 for PE, total reads for SE).
    """
    input:
        bam = "hisat2_aligned_processed/{name}_processed.bam"
    output:
        bdg = "scaled_bedgraph/{name}_processed.bdg"
    log: "logs/" + TIMESTAMP + "_{name}_make_scaled_bedgraph.log"
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
            echo '[{MODULE}] Generating RPM-scaled bedGraph for {{wildcards.name}}...' > {{log}}
            echo '[{MODULE}] Scaling factor: {scale_factor:.4f}' >> {{log}}

            bedtools genomecov \
                {pe} \
                -bga \
                -trackline \
                -trackopts 'name={{wildcards.name}} color=128,128,128' \
                -scale {scale_factor} \
                -ibam {{input.bam}} \
                > {{output.bdg}} \
                2>> {{log}}

            echo '[{MODULE}] ✓ bedGraph generation complete for {{wildcards.name}}' >> {{log}}
            '''
        )


# ================================================================================
# WORKFLOW COMPLETION
# ================================================================================

onsuccess:
    # Build completion summary
    layout = "paired-end" if PAIRED else "single-end"
    tn5_status = "Applied" if config['run_atacseq']['tn5_shift'] else "Not applied"
    exclusion_status = "Custom + mitochondria" if config['run_atacseq']['exclusion_list'] else "Mitochondria only"
    
    # Count total peaks across samples
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✓ ATAC-seq analysis complete!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Project: {PROJID}'
        echo '[{MODULE}] Samples: {len(NAMES)} ({layout})'
        echo '[{MODULE}] Tn5 shifting: {tn5_status}'
        echo '[{MODULE}] Exclusion regions: {exclusion_status}'
        echo ''
        echo '[{MODULE}] Completed analyses:'
        echo '[{MODULE}]   • Read preprocessing (fastp)'
        echo '[{MODULE}]   • Genome alignment (HISAT2, no splicing)'
        echo '[{MODULE}]   • Post-alignment filtering'
        if [ "{config['run_atacseq']['tn5_shift']}" = "True" ]; then
            echo '[{MODULE}]   • Tn5 transposase shifting'
        fi
        echo '[{MODULE}]   • Peak calling (MACS2)'
        echo '[{MODULE}]   • Peak quantification (featureCounts)'
        echo '[{MODULE}]   • Gene quantification (featureCounts)'
        echo '[{MODULE}]   • Coverage tracks (bedGraph, RPM-normalized)'
        echo ''
        
        # Count total peaks
        total_peaks=0
        for file in macs2_peaks/*_peaks.narrowPeak; do
            count=$(wc -l < "$file")
            total_peaks=$((total_peaks + count))
            sample=$(basename "$file" _processed_peaks.narrowPeak)
            echo "[{MODULE}]   • $sample: $count peaks"
        done
        echo "[{MODULE}]   • Total: $total_peaks peaks across all samples"
        echo ''
        
        echo '[{MODULE}] Output directories:'
        echo '[{MODULE}]   • fq_preprocessed/ (QC-filtered reads)'
        echo '[{MODULE}]   • hisat2_aligned/ (initial alignments)'
        echo '[{MODULE}]   • hisat2_aligned_processed/ (filtered BAM files)'
        echo '[{MODULE}]   • macs2_peaks/ (called peaks)'
        echo '[{MODULE}]   • featurecounts_reads_in_peaks/ (peak counts)'
        echo '[{MODULE}]   • featurecounts_reads_in_genes/ (gene counts)'
        echo '[{MODULE}]   • scaled_bedgraph/ (coverage tracks)'
        echo ''
        echo '[{MODULE}] Key output files:'
        echo '[{MODULE}]   • Peaks: macs2_peaks/*_processed_peaks.narrowPeak'
        echo '[{MODULE}]   • Peak counts: featurecounts_reads_in_peaks/*_reads_in_peaks.txt'
        echo '[{MODULE}]   • Gene counts: featurecounts_reads_in_genes/{PROJID}_reads_in_genes.txt'
        echo '[{MODULE}]   • Coverage: scaled_bedgraph/*_processed.bdg'
        echo ''
        echo '[{MODULE}] Next steps:'
        echo '[{MODULE}]   • Differential accessibility (DiffBind, DESeq2)'
        echo '[{MODULE}]   • Peak annotation and motif discovery (HOMER, MEME)'
        echo '[{MODULE}]   • TSS enrichment and QC metrics'
        echo '[{MODULE}]   • Genome browser visualization (IGV, UCSC)'
        echo ''
    """)

onerror:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✗ ATAC-seq analysis FAILED!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Common issues:'
        echo '[{MODULE}]   • Insufficient disk space for alignment files'
        echo '[{MODULE}]   • Memory errors during alignment or peak calling'
        echo '[{MODULE}]   • Missing or corrupted FASTQ input files'
        echo '[{MODULE}]   • Incorrect HISAT2 index path or corrupted index'
        echo '[{MODULE}]   • Missing chromosome sizes file'
        echo '[{MODULE}]   • Invalid exclusion list BED file format'
        echo '[{MODULE}]   • Incorrect MACS2 genome size parameter'
        echo '[{MODULE}]   • Low sequencing depth (<10M properly paired reads)'
        echo '[{MODULE}]   • GTF file format incompatible with featureCounts'
        echo '[{MODULE}]   • Required tools not in PATH'
        echo ''
        echo '[{MODULE}] Troubleshooting:'
        echo '[{MODULE}]   1. Check logs/ directory for detailed errors'
        echo '[{MODULE}]   2. Verify FASTQ files: {config["run_atacseq"]["indir"]}/*.f*q.gz'
        echo '[{MODULE}]   3. Confirm HISAT2 index: {config["run_atacseq"]["hisat2_ref"]}'
        echo '[{MODULE}]   4. Verify chromosome sizes: {config["run_atacseq"]["chr_sizes"]}'
        echo '[{MODULE}]   5. Check exclusion list: {config["run_atacseq"].get("exclusion_list", "None")}'
        echo '[{MODULE}]   6. Verify GTF file: {config["run_atacseq"]["gtf"]}'
        echo '[{MODULE}]   7. Confirm genome size: {config["run_atacseq"]["genome_size"]}'
        echo '[{MODULE}]   8. Check resources (CPU: {AVAIL_THREADS}, RAM, disk)'
        echo '[{MODULE}]   9. Verify sequencing depth is adequate'
        echo '[{MODULE}]   10. Ensure all tools installed: fastp, hisat2, samtools, MACS2, bedtools, deepTools, featureCounts'
        echo ''
        echo '[{MODULE}] Failed at timestamp: {TIMESTAMP}'
        echo ''
    """)
