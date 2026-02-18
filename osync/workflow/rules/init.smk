import os
import re
import glob

from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d")

# URLs for data retrieval
RIBORNA_LINK = 'https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz'
EXCLUSIONLIST_LINKS = [
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz'
]
KRAKEN2_STANDARD16_LINK = 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20221209.tar.gz'
KRAKEN2_PLUSPFP16_LINK = 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz'

# Set targets
RUN_RULES = []
INIT_HISAT2 = INIT_RIBORNA = INIT_EXCLUSIONLIST = INIT_CHRSIZESREF = INIT_KALLISTO = INIT_KRAKEN2 = None

if config["init"]["hisat2_ref"] and config["init"]["hisat2_index_base"]:
    INIT_HISAT2 = 'hisat2_index/' + config['init']['hisat2_index_base'] + '.1.ht2'
    RUN_RULES.append(INIT_HISAT2)

if config["init"]["rrna_ref"]:
    INIT_RIBORNA = 'rRNA_ref/rRNA_ref_complete'
    RUN_RULES.append(INIT_RIBORNA)

if config["init"]["exclusion_lists"]:
    INIT_EXCLUSIONLIST = 'exclusion_lists/exclusion_bed_files_complete'
    RUN_RULES.append(INIT_EXCLUSIONLIST)

if config["init"]["chr_sizes_ref"]:
    INIT_CHRSIZESREF = 'chromosome_sizes/' + os.path.basename(config["init"]["chr_sizes_ref"]) + '.sizes'
    RUN_RULES.append(INIT_CHRSIZESREF)

if config["init"]["kallisto_ref"] and config["init"]["kallisto_index"]:
    INIT_KALLISTO = 'kallisto_index/' + config['init']['kallisto_index']
    RUN_RULES.append(INIT_KALLISTO)

if config["init"]["kraken2_db_type"] and config["init"]["bracken_read_len"]:
    INIT_KRAKEN2 = (
        'kraken2_db/' + 
        config['init']['kraken2_db_type'] + 
        '/database' +
        str(config["init"]["bracken_read_len"]) +
        'mers.kmer_distrib'
    )
    RUN_RULES.append(INIT_KRAKEN2)

# Distribute allocated jobs/cores
AVAIL_THREADS = config['init']['jobs']
INIT_HISAT2_THREADS = INIT_RIBORNA_THREADS = INIT_EXCLUSIONLIST_THREADS = INIT_CHRSIZESREF_THREADS = INIT_KALLISTO_THREADS = INIT_KRAKEN2_THREADS = 0

if INIT_RIBORNA:
    INIT_RIBORNA_THREADS = 1
    AVAIL_THREADS -= 1

if INIT_EXCLUSIONLIST:
    INIT_EXCLUSIONLIST_THREADS = 1
    AVAIL_THREADS -= 1

if INIT_CHRSIZESREF:
    INIT_CHRSIZESREF_THREADS = 1
    AVAIL_THREADS -= 1

if INIT_KALLISTO:
    INIT_KALLISTO_THREADS = 1
    AVAIL_THREADS -= 1

if INIT_HISAT2:
    INIT_HISAT2_THREADS = round(AVAIL_THREADS/2) if INIT_KRAKEN2 else AVAIL_THREADS
    AVAIL_THREADS -= INIT_HISAT2_THREADS

if INIT_KRAKEN2:
    if config["init"]["kraken2_db_type"] in ['std16_prebuilt', 'plusPFP16_prebuilt']:
        INIT_KRAKEN2_THREADS = 1
        AVAIL_THREADS -= 1
    else:
        INIT_KRAKEN2_THREADS = AVAIL_THREADS
        AVAIL_THREADS = 0
    INIT_HISAT2_THREADS += AVAIL_THREADS


rule all:
    input:
        RUN_RULES


#
# Build HISAT2 index files
#
if INIT_HISAT2 and config['init']['hisat2_gtf']:

    # Include splice sites and exons
    rule hisat2_index:
        """
        Builds index for HISAT2.
        """
        input:
            ref = config['init']['hisat2_ref'],
            gtf = config['init']['hisat2_gtf']
        output:
            INIT_HISAT2
        params:
            ref = os.path.basename(config['init']['hisat2_ref']),
            gtf = os.path.basename(config['init']['hisat2_gtf']),
            base = config['init']['hisat2_index_base'],
            outdir = "hisat2_index"
        log: "logs/" + TIMESTAMP + "_hisat2_index.log"
        threads: INIT_HISAT2_THREADS
        # resources: cpus=4, mem_mb=4000, time_min=200
        shell:
            '''
            echo 'Extracting splicesites from '{params.gtf}'...' > {log}

            hisat2_extract_splice_sites.py \
                {input.gtf} \
                > {params.outdir}/{params.gtf}.ss.tsv \
                2>> {log}

            echo 'Finished extracting splicesites' >> {log}


            echo 'Extracting exons from '{params.gtf}'...' >> {log}

            hisat2_extract_exons.py \
                {input.gtf} \
                > {params.outdir}/{params.gtf}.exons.tsv \
                2>> {log}

            echo 'Finished extracting exons' >> {log}


            echo 'Building HISAT2 index using '{params.ref}'...' >> {log}

            hisat2-build \
                -p {threads} \
                --ss {params.outdir}/{params.gtf}.ss.tsv \
                --exon {params.outdir}/{params.gtf}.exons.tsv \
                {input.ref} {params.outdir}/{params.base} \
                &>> {log}

            echo 'Finished HISAT2 index' >> {log}
            '''

elif INIT_HISAT2:

    rule hisat2_index:
        """
        Builds index for HISAT2.
        """
        input:
            ref = config['init']['hisat2_ref']
        output:
            INIT_HISAT2
        params:
            ref = os.path.basename(config['init']['hisat2_ref']),
            base = config['init']['hisat2_index_base'],
            outdir = "hisat2_index"
        log: "logs/" + TIMESTAMP + "_hisat2_index.log"
        threads: INIT_HISAT2_THREADS
        # resources: cpus=4, mem_mb=4000, time_min=200
        shell:
            '''
            echo 'Building HISAT2 index using '{params.ref}'...' > {log}

            hisat2-build \
                -p {threads} \
                {input.ref} \
                {params.outdir}/{params.base} \
                &>> {log}

            echo 'Finished HISAT2 index' >> {log}
            '''


#
# Download rRNA reference files from SortMeRNA
#
if INIT_RIBORNA:

    rule download_rrna_ref:
            """
            Retrieves rRNA reference files.
            """
            output:
                complete = temp(INIT_RIBORNA)
            params:
                link = RIBORNA_LINK,
                filename = os.path.basename(RIBORNA_LINK),
                outdir = 'rRNA_ref'
            log: "logs/" + TIMESTAMP + "_download_rRNA_ref.log"
            threads: INIT_RIBORNA_THREADS
            # resources: cpus=4, mem_mb=4000, time_min=200
            shell:
                '''
                echo 'Downloading rRNA reference files...' > {log}

                curl -L {params.link} --remote-name --silent
                tar -zxf {params.filename} -C {params.outdir}
                rm {params.filename}

                echo 'Finished downloading rRNA reference' >> {log}

                touch {output.complete}
                '''


#
# Download reference exclusion list files
#
if INIT_EXCLUSIONLIST:

    rule download_exclusionlist_files:
            """
            Retrieves genomic region exclusion bed files.
            """
            output:
                complete = temp(INIT_EXCLUSIONLIST)
            params:
                links = EXCLUSIONLIST_LINKS,
                outdir = 'exclusion_lists'
            log: "logs/" + TIMESTAMP + "_download_exclusion_lists.log"
            threads: INIT_EXCLUSIONLIST_THREADS
            # resources: cpus=4, mem_mb=4000, time_min=200
            run:
                shell(" echo 'Downloading exclusion lists...' > {log} ")

                for link in params.links:
                    filename = os.path.basename(link)
                    new_filename = os.path.splitext(filename)[0]

                    shell(
                        '''
                        curl -L {link} --remote-name --silent
                        gunzip {filename}
                        sed -i 's/chr//g' {new_filename}
                        mv {new_filename} {params.outdir}/{new_filename}
                        '''
                    )
                
                shell(" echo 'Finished downloading exclusion lists' > {log} ")
                shell(" touch {output.complete} ")


#
# Retrieve chromosome sizes from reference genome
#
if INIT_CHRSIZESREF:

    rule get_chromosome_sizes:
            """
            Retrieves chromosome sizes.
            """
            input:
                config['init']['chr_sizes_ref']
            output:
                INIT_CHRSIZESREF
            params:
                ref = os.path.basename(config['init']['chr_sizes_ref']),
                outdir = "chromosome_sizes"
            log: "logs/" + TIMESTAMP + "_get_chromosome_sizes.log"
            threads: INIT_CHRSIZESREF_THREADS
            # resources: cpus=4, mem_mb=4000, time_min=200
            shell:
                '''
                echo 'Extracting chromosome sizes from reference file...' > {log}

                samtools faidx {input} &>> {log}
                awk '{{print $1"\t"$2}}' {input}.fai > {params.outdir}/{params.ref}.sizes
                rm {input}.fai

                echo 'Finished extracting chromosome sizes' >> {log}
                '''


#
# Builds kallisto index file
#
if INIT_KALLISTO:

    rule kallisto_index:
            """
            Builds index for kallisto.
            """
            input:
                ref = config['init']['kallisto_ref']
            output:
                idx = INIT_KALLISTO
            params:
                ref = os.path.basename(config['init']['kallisto_ref'])
            log: "logs/" + TIMESTAMP + "_kallisto_index.log"
            threads: INIT_KALLISTO_THREADS
            # resources: cpus=4, mem_mb=4000, time_min=200
            shell:
                '''
                echo 'Building kallisto index using '{params.ref}'...' > {log}

                kallisto index --index {output.idx} {input.ref} &>> {log}

                echo 'Finished kallisto index' >> {log}
                '''


#
# Builds database for kraken2/bracken analysis
#
if INIT_KRAKEN2:

    rule kraken2_db:
        """
        Builds kraken2/bracken database.
        """
        output:
            db = INIT_KRAKEN2
        params:
            db_type = config['init']['kraken2_db_type'],
            read_len = config['init']['bracken_read_len'],
            std16_link = KRAKEN2_STANDARD16_LINK,
            std16_filename = os.path.basename(KRAKEN2_STANDARD16_LINK),
            plusPFP16_link = KRAKEN2_PLUSPFP16_LINK,
            plusPFP16_filename = os.path.basename(KRAKEN2_PLUSPFP16_LINK),
            outdir = 'kraken2_db/' + config['init']['kraken2_db_type']
        log: "logs/" + TIMESTAMP + "_kraken2_db.log"
        threads: INIT_KRAKEN2_THREADS
        # resources: cpus=4, mem_mb=4000, time_min=200
        run:
            shell(" echo 'Building kraken2/bracken database for '{params.db_type}'...' > {log} ")

            if params['db_type'] == 'standard':
                shell(
                    '''
                    kraken2-build \
                        --threads {threads} \
                        --use-ftp \
                        --standard \
                        --db {params.outdir} \
                        &>> {log}
                    '''
                )
            elif params['db_type'] in ['greengenes','rdp','silva']:
                shell(
                    '''
                    kraken2-build \
                        --threads {threads} \
                        --use-ftp \
                        --special {params.db_type} \
                        --db {params.outdir} \
                        &>> {log}
                    '''
                )
            elif params['db_type'] == 'nt':
                shell(
                    '''
                    kraken2-build \
                        --use-ftp \
                        --download-taxonomy \
                        --db {params.outdir} \
                        &>> {log}

                    kraken2-build \
                        --use-ftp \
                        --download-library {params.db_type} \
                        --db {params.outdir} \
                        &>> {log}

                    kraken2-build \
                        --threads {threads} \
                        --build \
                        --fast-build \
                        --db {params.outdir} \
                        &>> {log}
                    '''
                )
            elif params['db_type'] == 'std16_prebuilt':
                shell(
                    '''
                    curl -L {params.std16_link} --remote-name --silent
                    tar -zxf {params.std16_filename} -C {params.outdir}
                    rm {params.std16_filename}

                    echo 'Finished kraken2/bracken database' >> {log}
                    '''
                )
            elif params['db_type'] == 'plusPFP16_prebuilt':
                shell(
                    '''
                    curl -L {params.plusPFP16_link} --remote-name --silent
                    tar -zxf {params.plusPFP16_filename} -C {params.outdir}
                    rm {params.plusPFP16_filename}

                    echo 'Finished kraken2/bracken database' >> {log}
                    '''
                )

            if params['db_type'] not in ['std16_prebuilt', 'plusPFP16_prebuilt']:
                shell(
                    '''
                    bracken-build \
                        -t {threads} \
                        -k 35 \
                        -l {params.read_len} \
                        -d {params.outdir} \
                        &>> {log}

                    echo 'Finished kraken2/bracken database' >> {log}
                    '''
                )


onsuccess:
    shell(" echo Initialization successful! ")
onerror:
    shell(" echo Initialization halted! ")
