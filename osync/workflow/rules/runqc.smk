import os
import re
import glob

from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d")

workflow_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
sys.path.append(os.path.dirname(workflow_dir))

from utils import check_strandedness, plot_biotype_composition

# Set variables
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

if config["runqc"]["hisat2_ref"] and config["runqc"]["gtf"]:
    HISAT2_THREADS = FASTP_THREADS

if config["runqc"]["kraken2_db"]:
    KRAKEN2_THREADS = AVAIL_THREADS
    BRACKEN_DB = glob.glob(os.path.join(config['runqc']['kraken2_db'], '*mers.kmer_distrib'))
    if len(BRACKEN_DB) == 1:
        BRACKEN_READ_LEN = int(re.search(r'database(\d+)mers.kmer_distrib', BRACKEN_DB[0]).group(1))
    else:
        BRACKEN_READ_LEN = 100

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


rule subset_fq:
    """
    Subsamples raw fastq files to desired read number/proportion.
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
        shell(
            '''
            echo 'Sampling '{wildcards.name}{wildcards.read}{EXT}' \
                (p={params.number})...' > {log}
            '''
        )

        if params['number'] < 1:
            shell(
                '''
                seqkit sample \
                    -j {threads} \
                    -p {params.number} \
                    -s {params.seed} \
                    -o {output.fq_subs} \
                    {input.fq} \
                    &>> {log}
                '''
            )

        elif params['number'] == 1:
            shell(" cp {input.fq} {output.fq_subs} ")

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

        shell(
            " echo 'Finished sampling '{wildcards.name}{wildcards.read}{EXT}'' >> {log} "
        )

rule fastqc_raw:
    """
    Runs fastqc on sampled raw fastq files.
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
    # resources: cpus=4, mem_mb=4000, time_min=200
    shell:
        '''
        echo 'Running QC on raw fastq files...' > {log}

        fastqc \
            --threads {threads} \
            --noextract \
            --outdir {params.outdir} \
            {input.fq_files} \
            &>> {log}

        echo 'Finished QC on raw fastq files' >> {log}
        '''

if PAIRED:
    rule preprocess_fq:
        """
        Performs preprocessing steps on raw fastq files.
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
        # resources: cpus=4, mem_mb=4000, time_min=200
        shell:
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
else:
    rule preprocess_fq:
        """
        Performs preprocessing steps on raw fastq files.
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
        # resources: cpus=4, mem_mb=4000, time_min=200
        shell:
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

rule fastqc_preprocessed:
    """
    Runs fastqc on sampled and preprocessed fastq files.
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
    # resources: cpus=4, mem_mb=4000, time_min=200
    shell:
        '''
        echo 'Running QC on preprocessed fastq files...' > {log}

        fastqc \
            --threads {threads} \
            --noextract \
            --outdir {params.outdir} \
            {input.fq_files} \
            &>> {log}

        echo 'Finished QC on preprocessed fastq files' >> {log}
        '''


#
# Check abundance of gene biotypes
#
if RUNQC_BIOTYPE_COMPOSITION:

    rule rseqc_bed:
        """
        Converts GTF to BED for RSeQC analysis.
        """
        input:
            gtf = config["runqc"]["gtf"]
        output:
            bed = "rseqc/" + GTF_BASENAME + "_transcripts.bed"
        params:
            outdir = "rseqc"
        log: "logs/" + TIMESTAMP + "_rseqc_bed.log"
        threads: 1
        # resources: cpus=4, mem_mb=4000, time_min=200
        shell:
            '''
            echo 'Converting GTF to BED for RSeQC analysis...' > {log}

            gtfToGenePred \
                {input.gtf} \
                {params.outdir}/transcript_genepred.tmp \
                &>> {log}
            
            genePredToBed \
                {params.outdir}/transcript_genepred.tmp \
                {output.bed} \
                &>> {log}

            rm {params.outdir}/transcript_genepred.tmp

            echo 'Finished GTF to BED conversion' >> {log}
            '''
    
    if PAIRED:
    
        rule rseqc_infer_experiment:
            """
            Evaluates strandedness of sequencing experiment with RSeQC. 
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
            # resources: cpus=20, mem_mb=40000, time_min=60
            shell:
                '''
                echo 'Inferring strandedness for '{wildcards.name}'...' > {log}

                # remove read unmapped (0x4), mate unmapped (0x8) and not primary alignment (0x100)
                hisat2 \
                    --threads {threads} \
                    -1 {input.fq1} \
                    -2 {input.fq2} \
                    -x {params.hisat2_ref} \
                    -S {params.outdir}/{wildcards.name}.sam \
                    &>> {log}

                samtools view \
                    -@ {threads} \
                    -bS \
                    -F 268 \
                    {params.outdir}/{wildcards.name}.sam \
                    > {params.outdir}/{wildcards.name}.bam \
                    2>> {log}
                
                samtools sort \
                    -@ {threads} \
                    {params.outdir}/{wildcards.name}.bam \
                    -o {params.outdir}/{wildcards.name}.bam \
                    &>> {log}
                
                samtools index \
                    -@ {threads} \
                    {params.outdir}/{wildcards.name}.bam \
                    &>> {log}

                infer_experiment.py \
                    -i {params.outdir}/{wildcards.name}.bam \
                    -r {input.bed} \
                    > {output.infer_experiment} \
                    2>> {log}

                rm {params.outdir}/{wildcards.name}.sam
                rm {params.outdir}/{wildcards.name}.bam
                rm {params.outdir}/{wildcards.name}.bam.bai

                echo 'Finished inferring strandedness for '{wildcards.name}'' >> {log}
                '''
        
        rule hisat2_align:
            """
            Uses hisat2 to align reads to a reference genome. 
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
            # resources: cpus=20, mem_mb=40000, time_min=60
            run:
                strandedness = check_strandedness(input['infer_experiment'], PAIRED)

                shell(" echo 'Inferred strandedness for '{wildcards.name}': {strandedness}' > {log} ")

                hisat2_strandedness = ''
                if strandedness == 'fwd':
                    hisat2_strandedness = '--rna-strandness FR'
                elif strandedness == 'rev':
                    hisat2_strandedness = '--rna-strandness RF'
                
                shell(
                    '''
                    echo 'HISAT2: aligning '{wildcards.name}'...' > {log}

                    hisat2 \
                        --threads {threads} \
                        {hisat2_strandedness} \
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
                )
    
    else:

        rule rseqc_infer_experiment:
            """
            Evaluates strandedness of sequencing experiment with RSeQC. 
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
            # resources: cpus=20, mem_mb=40000, time_min=60
            shell:
                '''
                echo 'Inferring strandedness for '{wildcards.name}'...' > {log}

                # remove read unmapped (0x4) and not primary alignment (0x100)
                hisat2 \
                    --threads {threads} \
                    -U {input.fq} \
                    -x {params.hisat2_ref} \
                    -S {params.outdir}/{wildcards.name}.sam \
                    &>> {log}

                samtools view \
                    -@ {threads} \
                    -bS \
                    -F 260 \
                    {params.outdir}/{wildcards.name}.sam \
                    > {params.outdir}/{wildcards.name}.bam \
                    2>> {log}
                
                samtools sort \
                    -@ {threads} \
                    {params.outdir}/{wildcards.name}.bam \
                    -o {params.outdir}/{wildcards.name}.bam \
                    &>> {log}
                
                samtools index \
                    -@ {threads} \
                    {params.outdir}/{wildcards.name}.bam \
                    &>> {log}

                infer_experiment.py \
                    -i {params.outdir}/{wildcards.name}.bam \
                    -r {input.bed} \
                    > {output.infer_experiment} \
                    2>> {log}

                rm {params.outdir}/{wildcards.name}.sam
                rm {params.outdir}/{wildcards.name}.bam
                rm {params.outdir}/{wildcards.name}.bam.bai

                echo 'Finished inferring strandedness for '{wildcards.name}'' >> {log}
                '''
        
        rule hisat2_align:
            """
            Uses hisat2 to align reads to a reference genome. 
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
            # resources: cpus=20, mem_mb=40000, time_min=60
            run:
                strandedness = check_strandedness(input['infer_experiment'], PAIRED)

                shell(" echo 'Inferred strandedness for '{wildcards.name}': {strandedness}' > {log} ")

                hisat2_strandedness = ''
                if strandedness == 'fwd':
                    hisat2_strandedness = '--rna-strandness F'
                elif strandedness == 'rev':
                    hisat2_strandedness = '--rna-strandness R'

                shell(
                    '''
                    echo 'HISAT2: aligning '{wildcards.name}'...' > {log}

                    hisat2 \
                        --threads {threads} \
                        {hisat2_strandedness} \
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
                )
    
    rule rseqc_read_distribution:
        """
        Checks distribution of mapped reads over genome features with RSeQC.
        """
        input:
            bam = "hisat2_aligned/{name}.bam",
            bed = "rseqc/" + GTF_BASENAME + "_transcripts.bed"
        output:
            read_distribution = "rseqc/{name}_rseqc_read_distribution.txt"
        log: "logs/" + TIMESTAMP + "_{name}_rseqc_read_distribution.log"
        threads: 1
        # resources: cpus=10, mem_mb=20000, time_min=300
        shell:
            '''
            echo 'Measuring read distribution for '{wildcards.name}'...' > {log}

            read_distribution.py \
                -i {input.bam} \
                -r {input.bed} \
                > {output.read_distribution} \
                2>> {log}

            echo 'Finished measuring read distribution for '{wildcards.name}'' >> {log}
            '''

    rule feature_counts:
        """
        Checks feature biotypes using featureCounts from subread.
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
        # resources: cpus=10, mem_mb=20000, time_min=300
        run:
            pe = '-p --countReadPairs -C -B' if PAIRED else ''

            strandedness = check_strandedness(input['infer_experiment'], PAIRED)

            shell(" echo 'Using inferred strandedness for '{NAMES[0]}': {strandedness}' > {log} ")

            featurecounts_strandedness = '-s 0'
            if strandedness == 'fwd':
                featurecounts_strandedness = '-s 1'
            elif strandedness == 'rev':
                featurecounts_strandedness = '-s 2'

            shell(
                '''
                echo 'Quantifying biotypes with featureCounts...' > {log}

                featureCounts \
                    -T {threads} \
                    {pe} \
                    {featurecounts_strandedness} \
                    -t exon \
                    -g {params.biotype} \
                    -a {input.gtf} \
                    -o {output.featurecounts} \
                    {input.bams} \
                    &>> {log}

                echo 'Finished biotype quantification with featureCounts' >> {log}
                '''
            )

            plot_biotype_composition(output.featurecounts, output.biotype_composition)


#
# Check abundance of different species
#
if RUNQC_SPECIES_COMPOSITION:

    if PAIRED:

        rule kraken2:
            """
            Performs taxonomic classification of reads with kraken2.
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
            # resources: cpus=4, mem_mb=60000, time_min=30
            shell:
                '''
                echo 'Running kraken2 for '{wildcards.name}'...' > {log}

                kraken2 \
                    --threads {threads} \
                    --paired \
                    --gzip-compressed \
                    --db {params.db} \
                    --report {output.report} \
                    --output {output.out} \
                    {input.fq1} \
                    {input.fq2} \
                    &>> {log}

                echo 'Finished kraken2 for '{wildcards.name}'' >> {log}
                '''

    else:

        rule kraken2:
            """
            Performs taxonomic classification of reads with kraken2. 
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
            # resources: cpus=4, mem_mb=60000, time_min=30
            shell:
                '''
                echo 'Running kraken2 for '{wildcards.name}'...' > {log}

                kraken2 \
                    --threads {threads} \
                    --gzip-compressed \
                    --db {params.db} \
                    --report {output.report} \
                    --output {output.out} \
                    {input.fq} \
                    &>> {log}

                echo 'Finished kraken2 for '{wildcards.name}'' >> {log}
                '''
    
    rule bracken:
        """
        Uses bracken to estimate species abundance from kraken2 results.
        """
        input:
            k2_report = "kraken2/{name}_kraken2_report.txt"
        output:
            out = "bracken/{name}_kraken2_bracken.txt",
            report = "bracken/{name}_kraken2_bracken_report.txt"
        params:
            db = config['runqc']['kraken2_db'],
            bracken_read_len = BRACKEN_READ_LEN
        log: "logs/" + TIMESTAMP + "_{name}_bracken.log"
        threads: KRAKEN2_THREADS
        # resources: cpus=4, mem_mb=2000, time_min=20
        shell:
            '''
            echo 'Running bracken on '{wildcards.name}'...' > {log}
            
            bracken \
                -r {params.bracken_read_len} \
                -l S \
                -t 5 \
                -d {params.db} \
                -i {input.k2_report} \
                -w {output.report} \
                -o {output.out} \
                &>> {log}
            
            echo 'Finished running bracken on '{wildcards.name}'' >> {log}
            '''

    rule krona:
        """
        Generates krona plots for species composition of samples.
        """
        input:
            reports = expand(["bracken/{name}_kraken2_bracken_report.txt"], name=NAMES)
        output:
            krona_html = "krona/" + PROJID + "_species_composition.html"
        params:
            outdir = "krona"
        log: "logs/" + TIMESTAMP + "_krona.log"
        threads: 1
        # resources: cpus=4, mem_mb=2000, time_min=20
        shell:
            '''
            echo 'Generating krona plots...' > {log}

            printf "%s\n" {NAMES} \
                | xargs -I{{}} \
                    cp bracken/{{}}_kraken2_bracken_report.txt bracken/{{}}_kraken2_bracken_report.tmp
            
            # fix for kreport2krona.py
            echo -e '0\t0\t0\tX\t0\t  X' | tee -a bracken/*_kraken2_bracken_report.tmp &> /dev/null

            printf "%s\n" {NAMES} \
                | xargs -I{{}} \
                    kreport2krona.py \
                        -r bracken/{{}}_kraken2_bracken_report.tmp \
                        -o {params.outdir}/{{}}.txt \
                        &>> {log}
            
            rm bracken/*_kraken2_bracken_report.tmp

            printf "%s\n" {NAMES} \
                | xargs -I{{}} \
                    sed -i.tmp -E 's/[a-z]__//g' {params.outdir}/{{}}.txt \
                    &>> {log}
            
            rm {params.outdir}/*.tmp

            ktImportText {params.outdir}/*.txt \
                -o {output.krona_html} \
                &>> {log}
            
            echo 'Finished krona plots' >> {log}
            '''


rule multiqc_report:
    """
    Generates combined quality control report.
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
    # resources: cpus=4, mem_mb=4000, time_min=200
    shell:
        '''
        echo 'Generating QC report...' > {log}

        multiqc \
            --force \
            -o {params.outdir} \
            {params.dirs} \
            &>> {log}

        echo 'Finished QC report' >> {log}
        '''


onsuccess:
    shell(" echo Execution of QC steps successful! ")
onerror:
    shell(" echo Execution of QC steps halted! ")
