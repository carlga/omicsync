import os
import re
import glob
import subprocess

from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d")

# Set variables
PROJID = os.path.basename(config['run_rnaseq']['indir'])
FILES = glob.glob(os.path.join(config['run_rnaseq']['indir'], '*.f*q.gz'))
EXT = re.search(r'[-\w]+(.f(ast)?q.gz)', FILES[0]).group(1)
try:
    NAMES = list(set([re.search(r'([-\w]+)(?=_[12])_[12].f(ast)?q.gz', x).group(1) for x in FILES]))
except AttributeError:
    NAMES = list(set([re.search(r'([-\w]+).f(ast)?q.gz', x).group(1) for x in FILES]))
PAIRED = False if len(NAMES) == len(FILES) else True
# READS = ['_1','_2'] if PAIRED else ['']

# Set cores and opts
AVAIL_THREADS = config['run_rnaseq']['jobs']
try:
    OPTS = config['run_rnaseq']['fastp_opt']
    FASTP_THREADS = int(re.search(r'(-w|--thread) (\d+)', OPTS).group(2))
    config['run_rnaseq']['fastp_opt'] = re.sub(r'(-w|--thread) (\d+) ','', OPTS)
except AttributeError:
    FASTP_THREADS = 3
FASTP_THREADS = AVAIL_THREADS if FASTP_THREADS > AVAIL_THREADS else FASTP_THREADS

try:
    OPTS = config['run_rnaseq']['hisat2_opt']
    HISAT2_THREADS = int(re.search(r'(-p|--thread) (\d+)', OPTS).group(2))
    config['run_rnaseq']['hisat2_opt'] = re.sub(r'(-p|--thread) (\d+) ', '', OPTS)
except AttributeError:
    HISAT2_THREADS = 3
HISAT2_THREADS = AVAIL_THREADS if HISAT2_THREADS > AVAIL_THREADS else HISAT2_THREADS

try:
    OPTS = config['run_rnaseq']['hisat2_opt']
    STRANDED = re.search(r'--rna-strandness(?= *(F|R)) *([FR][FR]?)', OPTS).group(2)
except AttributeError:
    STRANDED = False

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

if config["run_rnaseq"]["exon_quant"]:
    EXON_QUANT = True
else:
    EXON_QUANT = False

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


if PAIRED:

    rule preprocess_fq:
        """
        Performs preprocessing steps on raw fastq files.
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
            if params['removal_ref']:
                shell(
                    '''
                    echo 'Removing aligned reads from '{wildcards.name}'...' >> {log}

                    rm -rf {params.outdir}/sortmerna/{wildcards.name}*
                    rm -f {params.outdir}/sortmerna/out/{wildcards.name}*

                    sortmerna \
                        --threads {threads} \
                        --ref {params.removal_ref} \
                        --reads {output.fq1_preproc} \
                        --reads {output.fq2_preproc} \
                        --workdir {params.outdir}/sortmerna \
                        --kvdb {params.outdir}/sortmerna/{wildcards.name}_kvdb \
                        --readb {params.outdir}/sortmerna/{wildcards.name}_readb \
                        --fastx \
                        --aligned {params.outdir}/sortmerna/out/{wildcards.name}_removed \
                        --other {params.outdir}/{wildcards.name} \
                        --paired_in \
                        --out2 \
                        &>> {log}

                    mv {params.outdir}/sortmerna/out/{wildcards.name}_removed_fwd.fq.gz \
                        {params.outdir}/sortmerna/out/{wildcards.name}_removed_1.fastq.gz
                    mv {params.outdir}/sortmerna/out/{wildcards.name}_removed_rev.fq.gz \
                        {params.outdir}/sortmerna/out/{wildcards.name}_removed_2.fastq.gz

                    mv {params.outdir}/{wildcards.name}_fwd.fq.gz {output.fq1_preproc}
                    mv {params.outdir}/{wildcards.name}_rev.fq.gz {output.fq2_preproc}

                    mv {params.outdir}/sortmerna/out/{wildcards.name}_removed.log \
                        {params.outdir}/sortmerna/out/{wildcards.name}_sortmerna.log

                    echo 'Finished removing aligned reads from '{wildcards.name}'' >> {log}
                    '''
                )

else:

    rule preprocess_fq:
        """
        Performs preprocessing steps on raw fastq files.
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
            if params['removal_ref']:
                shell(
                    '''
                    echo 'Removing aligned reads from '{wildcards.name}'...' >> {log}
                    
                    rm -rf {params.outdir}/sortmerna/{wildcards.name}*
                    rm -f {params.outdir}/sortmerna/out/{wildcards.name}*

                    sortmerna \
                        --threads {threads} \
                        --ref {params.removal_ref} \
                        --reads {output.fq_preproc} \
                        --workdir {params.outdir}/sortmerna \
                        --kvdb {params.outdir}/sortmerna/{wildcards.name}_kvdb \
                        --readb {params.outdir}/sortmerna/{wildcards.name}_readb \
                        --fastx \
                        --aligned {params.outdir}/sortmerna/out/{wildcards.name}_removed \
                        --other {params.outdir}/{wildcards.name} \
                        &>> {log}

                    mv {params.outdir}/sortmerna/out/{wildcards.name}_removed.fq.gz \
                        {params.outdir}/sortmerna/out/{wildcards.name}_removed.fastq.gz

                    mv {params.outdir}/{wildcards.name}.fq.gz {output.fq_preproc}

                    mv {params.outdir}/sortmerna/out/{wildcards.name}_removed.log \
                        {params.outdir}/sortmerna/out/{wildcards.name}_sortmerna.log

                    echo 'Finished removing aligned reads from '{wildcards.name}'' >> {log}
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
            hisat2_ref = config['run_rnaseq']['hisat2_ref'],
            hisat2_opt = config['run_rnaseq']['hisat2_opt']
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

    if STRANDED:

        rule hisat2_split_bam:
            """
            Uses samtools to split bam into forward and reverse strand.
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
            # resources: cpus=20, mem_mb=40000, time_min=60
            run:
                if STRANDED == 'FR':
                    pair1 = "fwd"
                    pair2 = "rev"
                elif STRANDED == 'RF':
                    pair1 = "rev"
                    pair2 = "fwd"

                shell(
                    '''
                    echo 'Splitting bam for '{wildcards.name}'...' > {log}

                    # first strand
                    # include first in pair (0x40) and remove read reverse strand (0x10)
                    samtools view \
                        -@ {threads} \
                        -b \
                        -f 64 \
                        -F 16 \
                        {input.bam} \
                        > {params.outdir}/{wildcards.name}_{pair1}1.bam.tmp \
                        2>> {log}

                    # include second in pair (0x80) and read reverse strand (0x10)
                    samtools view \
                        -@ {threads} \
                        -b \
                        -f 144 \
                        {input.bam} \
                        > {params.outdir}/{wildcards.name}_{pair1}2.bam.tmp \
                        2>> {log}

                    samtools merge \
                        -@ {threads} \
                        -f \
                        {params.outdir}/{wildcards.name}_{pair1}1.bam.tmp \
                        {params.outdir}/{wildcards.name}_{pair1}2.bam.tmp \
                        -o {params.outdir}/{wildcards.name}_{pair1}.bam.tmp \
                        &>> {log}

                    samtools sort \
                        -@ {threads} \
                        {params.outdir}/{wildcards.name}_{pair1}.bam.tmp \
                        -o {params.outdir}/{wildcards.name}_{pair1}.bam \
                        &>> {log}

                    samtools index \
                        {params.outdir}/{wildcards.name}_{pair1}.bam \
                        &>> {log}


                    # second strand
                    # include second in pair (0x80) and remove read reverse strand (0x10)
                    samtools view \
                        -@ {threads} \
                        -b \
                        -f 128 \
                        -F 16 \
                        {input.bam} \
                        > {params.outdir}/{wildcards.name}_{pair2}1.bam.tmp \
                        2>> {log}

                    # include first in pair (0x40) and read reverse strand (0x10) 
                    samtools view \
                        -@ {threads} \
                        -b \
                        -f 80 \
                        {input.bam} \
                        > {params.outdir}/{wildcards.name}_{pair2}2.bam.tmp \
                        2>> {log}

                    samtools merge \
                        -@ {threads} \
                        -f \
                        {params.outdir}/{wildcards.name}_{pair2}1.bam.tmp \
                        {params.outdir}/{wildcards.name}_{pair2}2.bam.tmp \
                        -o {params.outdir}/{wildcards.name}_{pair2}.bam.tmp \
                        &>> {log}

                    samtools sort \
                        -@ {threads} \
                        {params.outdir}/{wildcards.name}_{pair2}.bam.tmp \
                        -o {params.outdir}/{wildcards.name}_{pair2}.bam \
                        &>> {log}
                    
                    samtools index \
                        {params.outdir}/{wildcards.name}_{pair2}.bam \
                        &>> {log}


                    rm {params.outdir}/{wildcards.name}*.tmp

                    echo 'Finished splitting bam for '{wildcards.name}'' >> {log}
                    '''
                )
    
    if RUN_RNASEQ_KALLISTO:

        rule kallisto_quant:
            """
            Quantifies transcript abundances with kallisto.
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
                '''
                echo 'kallisto: quantifying transcripts for '{wildcards.name}'...' > {log}

                kallisto quant \
                    --threads {threads} \
                    {params.kallisto_opt} \
                    --index {input.kallisto_ref} \
                    --output-dir {params.outdir}/{wildcards.name} \
                    {input.fq1} \
                    {input.fq2} \
                    &>> {log}
                
                mv {params.outdir}/{wildcards.name}/abundance.h5 {output.h5}
                mv {params.outdir}/{wildcards.name}/abundance.tsv {output.tsv}
                mv {params.outdir}/{wildcards.name}/run_info.json {output.json}

                rm -r {params.outdir}/{wildcards.name}

                echo 'kallisto: finished quantification for '{wildcards.name}'...' >> {log}
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
            hisat2_ref = config['run_rnaseq']['hisat2_ref'],
            hisat2_opt = config['run_rnaseq']['hisat2_opt']
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
    
    if STRANDED:

        rule hisat2_split_bam:
            """
            Uses samtools to split bam into forward and reverse strand.
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
            # resources: cpus=20, mem_mb=40000, time_min=60
            run:
                if STRANDED == 'F':
                    read1 = "fwd"
                    read2 = "rev"
                elif STRANDED == 'R':
                    read1 = "rev"
                    read2 = "fwd"

                shell(
                    '''
                    echo 'Splitting bam for '{wildcards.name}'...' > {log}

                    # first strand
                    # remove read reverse strand (0x10)
                    samtools view \
                        -@ {threads} \
                        -b \
                        -F 16 \
                        {input.bam} \
                        > {params.outdir}/{wildcards.name}_{read1}.bam \
                        2>> {log}

                    samtools index \
                        {params.outdir}/{wildcards.name}_{read1}.bam \
                        &>> {log}


                    # second strand
                    # include read reverse strand (0x10)
                    samtools view \
                        -@ {threads} \
                        -b \
                        -f 16 \
                        {input.bam} \
                        > {params.outdir}/{wildcards.name}_{read2}.bam \
                        2>> {log}
                    
                    samtools index \
                        {params.outdir}/{wildcards.name}_{read2}.bam \
                        &>> {log}

                    echo 'Finished splitting bam for '{wildcards.name}'' >> {log}
                    '''
                )

    if RUN_RNASEQ_KALLISTO:

        rule kallisto_quant:
            """
            Quantifies transcript abundances with kallisto.
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
                '''
                echo 'kallisto: quantifying transcripts for '{wildcards.name}'...' > {log}

                kallisto quant \
                    --threads {threads} \
                    {params.kallisto_opt} \
                    --single \
                    --fragment-length {params.frag_len} \
                    --sd {params.frag_len_sd} \
                    --index {input.kallisto_ref} \
                    --output-dir {params.outdir}/{wildcards.name} \
                    {input.fq} \
                    &>> {log}
                
                mv {params.outdir}/{wildcards.name}/abundance.h5 {output.h5}
                mv {params.outdir}/{wildcards.name}/abundance.tsv {output.tsv}
                mv {params.outdir}/{wildcards.name}/run_info.json {output.json}

                rm -r {params.outdir}/{wildcards.name}

                echo 'kallisto: finished quantification for '{wildcards.name}'...' >> {log}
                '''

rule feature_counts:
    """
    Gene expression quantification with featureCounts.
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
    # resources: cpus=10, mem_mb=20000, time_min=300
    run:
        pe = '-p --countReadPairs -C -B' if PAIRED else ''

        featurecounts_strandedness = '-s 0'
        if STRANDED in ['F','FR']:
            featurecounts_strandedness = '-s 1'
        elif STRANDED in ['R','RF']:
            featurecounts_strandedness = '-s 2'

        shell(
            '''
            echo 'Quantifying gene expression with featureCounts...' > {log}

            featureCounts \
                -T {threads} \
                {pe} \
                {featurecounts_strandedness} \
                -t exon \
                -g gene_id \
                -a {input.gtf} \
                -o {output.featurecounts} \
                {input.bams} \
                &>> {log}

            echo 'Finished gene expression quantification' >> {log}
            '''
        )

if EXON_QUANT:
    rule flatten_gtf:
        """
        Flatten GTF annotation for robust exon-level quantification.
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
        # resources: cpus=10, mem_mb=20000, time_min=300
        shell:
            '''
            echo 'Flattening GTF file...' > {log}

            flattenGTF -a {input.gtf} -o {output.saf_merge} -t exon -g gene_id \
                &>> {log}
            
            flattenGTF -a {input.gtf} -o {output.saf_chop} -t exon -g gene_id -C \
                &>> {log}

            echo 'Finished GTF file flattening' >> {log}
            '''
    
    rule feature_counts_exon:
        """
        Exon-level quantification with featureCounts.
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
                '''
                echo 'Quantifying exon expression with flattened GTF...' > {log}
                
                featureCounts \
                    -T {threads} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {input.saf_merge} \
                    -o {output.featurecounts_merge} \
                    --nonSplitOnly \
                    {input.bams} \
                    &>> {log}
                
                featureCounts \
                    -T {threads} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {input.saf_chop} \
                    -o {output.featurecounts_chop} \
                    --nonSplitOnly \
                    {input.bams} \
                    &>> {log}
                
                featureCounts \
                    -T {threads} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {input.saf_merge} \
                    -o {output.featurecounts_merge} \
                    -J \
                    {g} \
                    {input.bams} \
                    &>> {log}
                
                featureCounts \
                    -T {threads} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -F SAF \
                    -g GeneID \
                    -f \
                    -O \
                    -a {input.saf_chop} \
                    -o {output.featurecounts_chop} \
                    -J \
                    {g} \
                    {input.bams} \
                    &>> {log}
                
                echo 'Finished exon-level expression quantification' >> {log}
                '''
            )

rule make_bedgraph:
    """
    Generates genome coverage files for aligned reads. 
    """
    input:
        bam = "hisat2_aligned/{name}.bam"
    output:
        bdg = "hisat2_aligned_bedgraph/{name}.bdg"
    log: "logs/" + TIMESTAMP + "_{name}_make_bedgraph.log"
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
                -split \
                -trackline \
                -trackopts 'name={wildcards.name} color=128,128,128' \
                -scale {scale_factor} \
                -ibam {input.bam} \
                > {output.bdg} \
                2>> {log}

            echo 'Finished bedgraph for '{wildcards.name}'...' >> {log}
            '''
        )

if STRANDED:

    rule make_stranded_bedgraph:
        """
        Generates strand-specific genome coverage files for aligned reads. 
        """
        input:
            bam_fwd = "hisat2_aligned_split/{name}_fwd.bam",
            bam_rev = "hisat2_aligned_split/{name}_rev.bam"
        output:
            bdg_fwd = "hisat2_aligned_split_bedgraph/{name}_fwd.bdg",
            bdg_rev = "hisat2_aligned_split_bedgraph/{name}_rev.bdg"
        log: "logs/" + TIMESTAMP + "_{name}_make_stranded_bedgraph.log"
        threads: 1
        # resources: cpus=20, mem_mb=40000, time_min=60
        run:
            # calculate RPM scaling factor
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
                '''
                echo 'Generating strand-specific bedgraph for '{wildcards.name}'...' > {log}

                bedtools genomecov \
                    {pe} \
                    -bga \
                    -split \
                    -trackline \
                    -trackopts 'name={wildcards.name}_fwd color=255,30,30' \
                    -scale {scale_factor} \
                    -ibam {input.bam_fwd} \
                    > {output.bdg_fwd} \
                    2>> {log}
                
                bedtools genomecov \
                    {pe} \
                    -bga \
                    -split \
                    -trackline \
                    -trackopts 'name={wildcards.name}_rev color=30,30,255' \
                    -scale {scale_factor} \
                    -ibam {input.bam_rev} \
                    > {output.bdg_rev} \
                    2>> {log}

                echo 'Finished strand-specific bedgraph for '{wildcards.name}'...' >> {log}
                '''
            )


onsuccess:
    shell(" echo Execution of RNA-seq workflow successful! ")
onerror:
    shell(" echo Execution of RNA-seq workflow halted! ")
