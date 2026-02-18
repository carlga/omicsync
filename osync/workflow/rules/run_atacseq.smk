import os
import re
import glob
import subprocess

from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d")

# Set variables
PROJID = os.path.basename(config['run_atacseq']['indir'])
FILES = glob.glob(os.path.join(config['run_atacseq']['indir'], '*.f*q.gz'))
EXT = re.search(r'[-\w]+(.f(ast)?q.gz)', FILES[0]).group(1)
try:
    NAMES = list(set([re.search(r'([-\w]+)(?=_[12])_[12].f(ast)?q.gz', x).group(1) for x in FILES]))
except AttributeError:
    NAMES = list(set([re.search(r'([-\w]+).f(ast)?q.gz', x).group(1) for x in FILES]))
PAIRED = False if len(NAMES) == len(FILES) else True
# READS = ['_1','_2'] if PAIRED else ['']

# Set cores and opts
AVAIL_THREADS = config['run_atacseq']['jobs']
try:
    OPTS = config['run_atacseq']['fastp_opt']
    FASTP_THREADS = int(re.search(r'(-w|--thread) (\d+)', OPTS).group(2))
    config['run_atacseq']['fastp_opt'] = re.sub(r'(-w|--thread) (\d+) ','', OPTS)
except AttributeError:
    FASTP_THREADS = 3
FASTP_THREADS = AVAIL_THREADS if FASTP_THREADS > AVAIL_THREADS else FASTP_THREADS

try:
    OPTS = config['run_atacseq']['hisat2_opt']
    HISAT2_THREADS = int(re.search(r'(-p|--thread) (\d+)', OPTS).group(2))
    config['run_atacseq']['hisat2_opt'] = re.sub(r'(-p|--thread) (\d+) ', '', OPTS)
except AttributeError:
    HISAT2_THREADS = 3
HISAT2_THREADS = AVAIL_THREADS if HISAT2_THREADS > AVAIL_THREADS else HISAT2_THREADS

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


if PAIRED:

    rule preprocess_fq:
        """
        Performs preprocessing steps on raw fastq files.
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
        # resources: cpus=4, mem_mb=4000, time_min=200
        run:
            shell(
                '''
                echo 'Processing raw reads from '{wildcards.name}'...' > {log}

                fastp \
                    --thread {threads} \
                    {params.fastp_opt} \
                    -i {input.fq1} \
                    -o {output.fq1_preproc} \
                    -I {input.fq2} \
                    -O {output.fq2_preproc} \
                    -j {output.json} \
                    -h {output.html} \
                    &>> {log}

                echo 'Finished processing reads from '{wildcards.name}'' >> {log}
                '''
            )

else:

    rule preprocess_fq:
        """
        Performs preprocessing steps on raw fastq files.
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
        # resources: cpus=4, mem_mb=4000, time_min=200
        run:
            shell(
                '''
                echo 'Processing raw reads from '{wildcards.name}'...' > {log}

                fastp \
                    --thread {threads} \
                    {params.fastp_opt} \
                    -i {input.fq} \
                    -o {output.fq_preproc} \
                    -j {output.json} \
                    -h {output.html} \
                    &>> {log}

                echo 'Finished processing reads from '{wildcards.name}'' >> {log}
                '''
            )

if PAIRED:

    rule hisat2_align:
        """
        Uses hisat2 to align reads to a reference genome. 
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
        # resources: cpus=20, mem_mb=40000, time_min=60
        shell:
            '''
            echo 'HISAT2: aligning '{wildcards.name}'...' > {log}

            hisat2 \
                --threads {threads} \
                {params.hisat2_opt} \
                -1 {input.fq1} \
                -2 {input.fq2} \
                -x {params.hisat2_ref} \
                --no-spliced-alignment \
                --no-temp-splicesite \
                --summary-file {output.summary} \
                -S {output.sam} \
                &>> {log}

            echo 'HISAT2: finished aligning '{wildcards.name}'' >> {log}

            # remove read unmapped (0x4), mate unmapped (0x8) and not primary alignment (0x100)
            samtools view \
                -@ {threads} \
                -bS \
                -F 268 \
                {output.sam} \
                > {output.bam} \
                2>> {log}
            
            samtools sort \
                -@ {threads} \
                {output.bam} \
                -o {output.bam} \
                &>> {log}
            
            samtools index \
                -@ {threads} \
                {output.bam} \
                &>> {log}

            samtools flagstat \
                -@ {threads} \
                {output.bam} \
                > {output.stats} \
                2>> {log}

            echo 'HISAT2: generated bam and stats for '{wildcards.name}'' >> {log}
            '''

else:

    rule hisat2_align:
        """
        Uses hisat2 to align reads to a reference genome. 
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
        # resources: cpus=20, mem_mb=40000, time_min=60
        shell:
            '''
            echo 'HISAT2: aligning '{wildcards.name}'...' > {log}

            hisat2 \
                --threads {threads} \
                {params.hisat2_opt} \
                -U {input.fq} \
                -x {params.hisat2_ref} \
                --no-spliced-alignment \
                --no-temp-splicesite \
                --summary-file {output.summary} \
                -S {output.sam} \
                &>> {log}

            echo 'HISAT2: finished aligning '{wildcards.name}'' >> {log}

            # remove read unmapped (0x4) and not primary alignment (0x100)
            samtools view \
                -@ {threads} \
                -bS \
                -F 260 \
                {output.sam} \
                > {output.bam} \
                2>> {log}
            
            samtools sort \
                -@ {threads} \
                {output.bam} \
                -o {output.bam} \
                &>> {log}
            
            samtools index \
                -@ {threads} \
                {output.bam} \
                &>> {log}

            samtools flagstat \
                -@ {threads} \
                {output.bam} \
                > {output.stats} \
                2>> {log}

            echo 'HISAT2: generated bam and stats for '{wildcards.name}'' >> {log}
            '''

rule generate_inclusion_list:
    """
    Generates list of genomic regions to include.
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
    # resources: cpus=10, mem_mb=20000, time_min=300
    run:
        if input.exclusion_list:
            shell(
                '''
                echo 'Generating inclusion list...' > {log}

                awk '$1 ~ /(MT|mt|chrMT|chrmt|chrM|chrm)/ {{print $1"\t0\t"$2}}' {input.chr_sizes} \
                    > {params.outdir}/exclusion_list.tmp \
                    2>> {log}
                
                awk '{{print $1,$2,$3}}' OFS='\t' {input.exclusion_list} \
                    >> {params.outdir}/exclusion_list.tmp \
                    2>> {log}

                sortBed -g {input.chr_sizes} -i {params.outdir}/exclusion_list.tmp \
                    > {params.outdir}/exclusion_list_sorted.tmp \
                    2>> {log}
                
                complementBed -g {input.chr_sizes} -i {params.outdir}/exclusion_list_sorted.tmp \
                    > {params.outdir}/inclusion_list.bed \
                    2>> {log}
                
                rm {params.outdir}/*.tmp

                echo 'Finished generating inclusion list' >> {log}
                '''
            )
        else:
            shell(
                '''
                echo 'Generating inclusion list...' > {log}

                awk '$1 ~ /(MT|mt|chrMT|chrmt|chrM|chrm)/ {{print $1"\t0\t"$2}}' {input.chr_sizes} \
                    > {params.outdir}/exclusion_list.tmp \
                    2>> {log}

                sortBed -g {input.chr_sizes} -i {params.outdir}/exclusion_list.tmp \
                    > {params.outdir}/exclusion_list_sorted.tmp \
                    2>> {log}
                
                complementBed -g {input.chr_sizes} -i {params.outdir}/exclusion_list_sorted.tmp \
                    > {params.outdir}/inclusion_list.bed \
                    2>> {log}
                
                rm {params.outdir}/*.tmp

                echo 'Finished generating inclusion list' >> {log}
                '''
            )

rule hisat2_align_process:
    """
    Filters and shifts aligned reads.
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
    # resources: cpus=10, mem_mb=20000, time_min=300
    run:
        shell(" echo 'Processing aligned reads...' > {log} ")

        shell(
            '''
            samtools view \
                -@ {threads} \
                -b \
                -L {input.inclusion_list} \
                {input.bam} \
                > {output.processed_bam} \
                2>> {log}
            
            samtools sort \
                -@ {threads} \
                {output.processed_bam} \
                -o {output.processed_bam} \
                &>> {log}
            
            samtools index \
                -@ {threads} \
                {output.processed_bam} \
                &>> {log}
            '''
        )

        if config['run_atacseq']['tn5_shift']:
            shell(
                '''
                alignmentSieve \
                    -p {threads} \
                    --ATACshift \
                    -b {output.processed_bam} \
                    -o {output.processed_bam} \
                    &>> {log}

                samtools sort \
                    -@ {threads} \
                    {output.processed_bam} \
                    -o {output.processed_bam} \
                    &>> {log}
                
                samtools index \
                    -@ {threads} \
                    {output.processed_bam} \
                    &>> {log}
                '''
            )
        shell(" echo 'Finished processing aligned reads' > {log} ")

rule macs2_call_peaks:
    """
    Calls peaks with MACS2.
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
    # resources: cpus=10, mem_mb=20000, time_min=300
    run:
        fmt = 'BAMPE' if PAIRED else 'BAM'

        shell(
            '''
            echo 'Calling peaks with MACS2...' > {log}

            macs2 callpeak \
                --bdg \
                -f {fmt} \
                -g {params.genome_size} \
                -t {input.bam} \
                -n {wildcards.name}_processed \
                --outdir {params.outdir} \
                &>> {log}

            echo 'Finished calling peaks' >> {log}
            '''
        )

rule reads_in_peaks:
    """
    Quantifies reads/fragments falling into called peak regions.
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
    # resources: cpus=10, mem_mb=20000, time_min=300
    run:
        pe = '-p --countReadPairs -C -B' if PAIRED else ''

        shell(
            '''
            echo 'Quantifying reads/fragments in called peaks with featureCounts...' > {log}

            echo -e 'GeneID\tChr\tStart\tEnd\tStrand' > {output.saf} 2>> {log}
            awk '{{print $4,$1,$2+1,$3,"."}}' OFS='\t' {input.peaks_bed} \
                >> {output.saf} \
                2>> {log}

            featureCounts \
                -T {threads} \
                {pe} \
                -F SAF \
                -a {output.saf} \
                -o {output.reads_in_peaks} \
                {input.bam} \
                &>> {log}

            echo 'Finished reads/fragment quantification' >> {log}
            '''
        )

rule reads_in_genes:
    """
    Quantifies reads/fragments falling into annotated genes.
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
    # resources: cpus=10, mem_mb=20000, time_min=300
    run:
        pe = '-p --countReadPairs -C -B' if PAIRED else ''

        shell(
            '''
            echo 'Quantifying reads/fragments in annotated genes with featureCounts...' > {log}

            featureCounts \
                -T {threads} \
                {pe} \
                -t gene \
                -g gene_id \
                -a {input.gtf} \
                -o {output.reads_in_genes} \
                {input.bams} \
                &>> {log}

            echo 'Finished reads/fragment quantification' >> {log}
            '''
        )

rule make_scaled_bedgraph:
    """
    Generates genome coverage files for aligned reads. 
    """
    input:
        bam = "hisat2_aligned_processed/{name}_processed.bam"
    output:
        bdg = "scaled_bedgraph/{name}_processed.bdg"
    log: "logs/" + TIMESTAMP + "_{name}_make_scaled_bedgraph.log"
    threads: 1
    # resources: cpus=20, mem_mb=40000, time_min=60
    run:
        # calculate RPM scaling factor
        reads = subprocess.check_output(
            f"samtools view -c {input.bam}", 
            shell=True, 
            text=True
        ).strip()
        scale_factor = 1e6/(int(reads)/2) if PAIRED else 1e6/int(reads)

        pe = '-pc' if PAIRED else ''

        shell(
            '''
            echo 'Generating bedgraph for '{wildcards.name}'...' > {log}

            bedtools genomecov \
                {pe} \
                -bga \
                -trackline \
                -trackopts 'name={wildcards.name} color=128,128,128' \
                -scale {scale_factor} \
                -ibam {input.bam} \
                > {output.bdg} \
                2>> {log}

            echo 'Finished bedgraph for '{wildcards.name}'...' >> {log}
            '''
        )


onsuccess:
    shell(" echo Execution of ATAC-seq workflow successful! ")
onerror:
    shell(" echo Execution of ATAC-seq workflow halted! ")
