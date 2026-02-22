"""
================================================================================
INIT MODULE - Reference Database Initialization Workflow
================================================================================
Initializes reference databases and index files for downstream analysis workflows.

Supported components:
- HISAT2 index building (with/without GTF annotation for splice-aware alignment)
- Kallisto index building (for transcript quantification)
- Kraken2/Bracken databases (taxonomic classification)
- rRNA reference databases (for contamination removal with SortMeRNA)
- Genomic exclusion lists (blacklisted regions for ATAC-seq/ChIP-seq)
- Chromosome size files (for genome coverage analysis)

This workflow handles both building from source and downloading pre-built databases.
All pre-built RefSeq databases include Bracken for 50, 75, 100, 150, 200, 250, 300-mers.
All pre-built 16S databases include Bracken for 100, 150, 200-mers.
================================================================================
"""

import os
import sys
from pathlib import Path

# Setup workflow environment
workflow_dir = Path(workflow.snakefile).parent.absolute()
sys.path.insert(0, str(workflow_dir.parent))

# Import utilities
from utils.workflow_utils import get_timestamp, allocate_threads

# Module identifier
MODULE = 'init'
TIMESTAMP = get_timestamp()


# ================================================================================
# DATA SOURCE URLS
# ================================================================================

# SortMeRNA rRNA reference database (v4.3.4)
RIBORNA_LINK = 'https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz'

# ENCODE Blacklist genomic exclusion regions (v2)
EXCLUSIONLIST_LINKS = [
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz'
]

# Pre-built Kraken2 databases (updated October 2025)
# Format: {db_name: (url, archive_size_GB, index_size_GB, description)}
KRAKEN2_PREBUILT_DBS = {
    # Viral databases
    'viral': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz',
        0.5, 0.6,
        'Refseq viral'
    ),
    
    # MinusB - No bacteria
    'minusb': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20251015.tar.gz',
        7.8, 11.0,
        'Refseq archaea, viral, plasmid, human, UniVec_Core (no bacteria)'
    ),
    
    # Standard databases
    'standard': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20251015.tar.gz',
        72.5, 94.0,
        'Refseq archaea, bacteria, viral, plasmid, human, UniVec_Core'
    ),
    'standard_8': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz',
        5.5, 7.5,
        'Standard with DB capped at 8 GB'
    ),
    'standard_16': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20251015.tar.gz',
        11.2, 14.9,
        'Standard with DB capped at 16 GB'
    ),
    
    # PlusPF - Standard + Protozoa & Fungi
    'pluspf': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20251015.tar.gz',
        77.5, 100.6,
        'Standard plus Refseq protozoa & fungi'
    ),
    'pluspf_8': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20251015.tar.gz',
        5.5, 7.5,
        'PlusPF with DB capped at 8 GB'
    ),
    'pluspf_16': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20251015.tar.gz',
        11.2, 14.9,
        'PlusPF with DB capped at 16 GB'
    ),
    
    # PlusPFP - Standard + Protozoa, Fungi & Plant
    'pluspfp': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20251015.tar.gz',
        158.8, 214.5,
        'Standard plus Refseq protozoa, fungi & plant'
    ),
    'pluspfp_8': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08_GB_20251015.tar.gz',
        5.2, 7.5,
        'PlusPFP with DB capped at 8 GB'
    ),
    'pluspfp_16': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16_GB_20251015.tar.gz',
        10.5, 14.9,
        'PlusPFP with DB capped at 16 GB'
    ),
    
    # 16S rRNA databases (includes Bracken for 100/150/200mers)
    '16s_greengenes_13.5': (
        'https://genome-idx.s3.amazonaws.com/kraken/16S_Greengenes13.5_20200326.tgz',
        0.073, 0.073,
        'Greengenes 13.5 16S rRNA (with Bracken 100/150/200mers)'
    ),
    '16s_rdp_11.5': (
        'https://genome-idx.s3.amazonaws.com/kraken/16S_RDP11.5_20200326.tgz',
        0.168, 0.168,
        'RDP 11.5 16S rRNA (with Bracken 100/150/200mers)'
    ),
    '16s_silva_132': (
        'https://genome-idx.s3.amazonaws.com/kraken/16S_Silva132_20200326.tgz',
        0.117, 0.117,
        'Silva 132 16S rRNA (with Bracken 100/150/200mers)'
    ),
    '16s_silva_138': (
        'https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz',
        0.112, 0.112,
        'Silva 138 16S rRNA (with Bracken 100/150/200mers)'
    ),
    
    # PrackenDB - Single reference per species
    'pracken': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_pracken_20251007.tar.gz',
        351.6, 351.6,
        'PrackenDB: Single reference per species (bacteria, archaea, protists, fungi, human, viral)'
    ),
    
    # Specialized databases
    'eupathdb48': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20230407.tar.gz',
        8.4, 11.0,
        'Eukaryotic pathogen genomes with contaminants removed'
    ),
    'core_nt': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz',
        245.8, 316.2,
        'Very large collection: GenBank, RefSeq, TPA and PDB'
    ),
    'gtdb_v226': (
        'https://genome-idx.s3.amazonaws.com/kraken/k2_gtdb_genome_reps_20250609.tar.gz',
        500.9, 644.0,
        'GTDB v226 bacterial and archaeal genomes'
    )
}


# ================================================================================
# CONFIGURATION AND SETUP
# ================================================================================

# Determine targets based on configuration
RUN_RULES = []
INIT_HISAT2 = INIT_RIBORNA = INIT_EXCLUSIONLIST = None
INIT_CHRSIZESREF = INIT_KALLISTO = INIT_KRAKEN2 = None

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

# For Kraken2: only require bracken_read_len for custom builds
if config["init"]["kraken2_db_type"]:
    db_type = config['init']['kraken2_db_type']
    
    # Pre-built databases: create a marker file (don't need specific kmer file)
    if db_type in KRAKEN2_PREBUILT_DBS:
        INIT_KRAKEN2 = f'kraken2_db/{db_type}/database_download_complete'
    # Custom builds: create the specific Bracken kmer distribution file
    else:
        if config["init"]["bracken_read_len"]:
            INIT_KRAKEN2 = (
                f'kraken2_db/{db_type}/' + 
                f'database{config["init"]["bracken_read_len"]}mers.kmer_distrib'
            )
        else:
            # Default to 100 for custom builds if not specified
            INIT_KRAKEN2 = f'kraken2_db/{db_type}/database100mers.kmer_distrib'
    
    RUN_RULES.append(INIT_KRAKEN2)


# ================================================================================
# WILDCARD CONSTRAINTS
# ================================================================================

wildcard_constraints:
    filename=r"[^/]+\.gz"


# ================================================================================
# THREAD ALLOCATION
# ================================================================================

AVAIL_THREADS = config['init']['jobs']
task_priorities = {}

# Priority 1.0: Simple download/extraction tasks (single-threaded)
if INIT_RIBORNA:
    task_priorities['riborna'] = 1
if INIT_EXCLUSIONLIST:
    task_priorities['exclusionlist'] = 1
if INIT_CHRSIZESREF:
    task_priorities['chrsizes'] = 1
if INIT_KALLISTO:
    task_priorities['kallisto'] = 1

# Priority 0.4: Compute-intensive building tasks (multi-threaded)
if INIT_HISAT2:
    task_priorities['hisat2'] = 0.4
if INIT_KRAKEN2:
    # Pre-built databases: just download (single-threaded)
    if config["init"]["kraken2_db_type"] in KRAKEN2_PREBUILT_DBS:
        task_priorities['kraken2'] = 1
    # Custom builds: multi-threaded
    else:
        task_priorities['kraken2'] = 0.4

# Allocate threads based on priorities
thread_allocation = allocate_threads(AVAIL_THREADS, task_priorities)
INIT_HISAT2_THREADS = thread_allocation.get('hisat2', 0)
INIT_RIBORNA_THREADS = thread_allocation.get('riborna', 0)
INIT_EXCLUSIONLIST_THREADS = thread_allocation.get('exclusionlist', 0)
INIT_CHRSIZESREF_THREADS = thread_allocation.get('chrsizes', 0)
INIT_KALLISTO_THREADS = thread_allocation.get('kallisto', 0)
INIT_KRAKEN2_THREADS = thread_allocation.get('kraken2', 0)


# ================================================================================
# WORKFLOW INITIALIZATION
# ================================================================================

onstart:
    # Build enabled tasks list
    enabled_tasks = []
    
    if INIT_HISAT2:
        gtf_status = "with GTF" if config["init"].get("hisat2_gtf") else "no GTF"
        enabled_tasks.append(f"HISAT2 index ({gtf_status})")
    
    if INIT_KALLISTO:
        enabled_tasks.append("Kallisto index")
    
    if INIT_KRAKEN2:
        db_type = config['init']['kraken2_db_type']
        if db_type in KRAKEN2_PREBUILT_DBS:
            _, archive_size, index_size, desc = KRAKEN2_PREBUILT_DBS[db_type]
            if db_type.startswith('16s_'):
                enabled_tasks.append(f"Kraken2/Bracken: {db_type} ({archive_size} GB, Bracken 100/150/200mers)")
            else:
                enabled_tasks.append(f"Kraken2/Bracken: {db_type} ({index_size} GB, Bracken 50-300mers)")
        else:
            read_len = config.get('init', {}).get('bracken_read_len', 100)
            enabled_tasks.append(f"Kraken2/Bracken: {db_type} (custom build, Bracken {read_len}mers)")
    
    if INIT_RIBORNA:
        enabled_tasks.append("rRNA references (SortMeRNA v4.3.4)")
    
    if INIT_EXCLUSIONLIST:
        enabled_tasks.append(f"Exclusion lists ({len(EXCLUSIONLIST_LINKS)} genomes)")
    
    if INIT_CHRSIZESREF:
        ref_name = os.path.basename(config["init"]["chr_sizes_ref"])
        enabled_tasks.append(f"Chromosome sizes ({ref_name})")
    
    tasks_str = "\n".join([f"[{MODULE}]   • {task}" for task in enabled_tasks])
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Initialization Workflow'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        if [ {len(enabled_tasks)} -eq 0 ]; then
            echo '[{MODULE}] ⚠ No initialization tasks configured!'
        else
            echo '[{MODULE}] Tasks: {len(enabled_tasks)}'
            echo ''
            echo '[{MODULE}] Enabled components:'
            echo '{tasks_str}'
        fi
        echo ''
        echo '[{MODULE}] Threads: {AVAIL_THREADS}'
        echo '[{MODULE}] Timestamp: {TIMESTAMP}'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo ''
    """)


# ================================================================================
# TARGET RULES
# ================================================================================

rule all:
    input: RUN_RULES


# ================================================================================
# HISAT2 INDEX BUILDING
# ================================================================================

if INIT_HISAT2 and config['init']['hisat2_gtf']:
    rule hisat2_index:
        """
        Builds HISAT2 index with splice-aware alignment support.
        Extracts splice sites and exons from GTF annotation for improved
        alignment accuracy in RNA-seq workflows. Suitable for transcriptome analysis.
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
        log: f"logs/{TIMESTAMP}_hisat2_index.log"
        threads: INIT_HISAT2_THREADS
        shell:
            f'''
            echo '[{MODULE}] Extracting splice sites from {{params.gtf}}...' | tee {{log}}
            hisat2_extract_splice_sites.py {{input.gtf}} > {{params.outdir}}/{{params.gtf}}.ss.tsv 2>> {{log}}

            echo '[{MODULE}] Extracting exons from {{params.gtf}}...' | tee -a {{log}}
            hisat2_extract_exons.py {{input.gtf}} > {{params.outdir}}/{{params.gtf}}.exons.tsv 2>> {{log}}

            echo '[{MODULE}] Building HISAT2 index with splice-aware annotation...' | tee -a {{log}}
            hisat2-build -p {{threads}} \
                --ss {{params.outdir}}/{{params.gtf}}.ss.tsv \
                --exon {{params.outdir}}/{{params.gtf}}.exons.tsv \
                {{input.ref}} {{params.outdir}}/{{params.base}} \
                2>> {{log}}

            echo '[{MODULE}] ✓ HISAT2 index complete (with GTF annotation)' | tee -a {{log}}
            '''

elif INIT_HISAT2:
    rule hisat2_index:
        """
        Builds basic HISAT2 index from reference genome without annotation.
        Suitable for DNA alignment (ATAC-seq, ChIP-seq) or when GTF is unavailable.
        """
        input:
            ref = config['init']['hisat2_ref']
        output:
            INIT_HISAT2
        params:
            base = config['init']['hisat2_index_base'],
            outdir = "hisat2_index"
        log: f"logs/{TIMESTAMP}_hisat2_index.log"
        threads: INIT_HISAT2_THREADS
        shell:
            f'''
            echo '[{MODULE}] Building HISAT2 index from reference genome...' | tee {{log}}
            hisat2-build -p {{threads}} {{input.ref}} {{params.outdir}}/{{params.base}} 2>> {{log}}
            echo '[{MODULE}] ✓ HISAT2 index complete (no GTF annotation)' | tee -a {{log}}
            '''


# ================================================================================
# rRNA REFERENCE DOWNLOAD
# ================================================================================

if INIT_RIBORNA:
    rule download_rrna_ref:
        """
        Downloads rRNA reference database from SortMeRNA v4.3.4.
        Used for removing ribosomal RNA contamination from RNA-seq data.
        Contains curated rRNA sequences from multiple organisms.
        """
        output:
            complete = temp(INIT_RIBORNA)
        params:
            link = RIBORNA_LINK,
            filename = os.path.basename(RIBORNA_LINK),
            outdir = 'rRNA_ref'
        log: f"logs/{TIMESTAMP}_download_rRNA_ref.log"
        threads: INIT_RIBORNA_THREADS
        shell:
            f'''
            echo '[{MODULE}] Downloading rRNA reference database (SortMeRNA v4.3.4)...' | tee {{log}}
            curl -L {{params.link}} --remote-name --silent
            
            echo '[{MODULE}] Extracting rRNA database...' | tee -a {{log}}
            tar -zxf {{params.filename}} -C {{params.outdir}}
            rm {{params.filename}}
            
            echo '[{MODULE}] ✓ rRNA references downloaded and extracted' | tee -a {{log}}
            touch {{output.complete}}
            '''


# ================================================================================
# EXCLUSION LIST DOWNLOAD
# ================================================================================

if INIT_EXCLUSIONLIST:
    rule download_exclusionlist_files:
        """
        Downloads ENCODE blacklist BED files for multiple genome assemblies.
        These regions are known to produce artifacts in high-throughput sequencing
        (ATAC-seq, ChIP-seq, etc.) and should be excluded from peak calling.
        
        Supported genomes: ce10, ce11, dm3, dm6, hg19, hg38, mm10
        Removes 'chr' prefix for compatibility with Ensembl-style chromosome names.
        """
        output:
            complete = temp(INIT_EXCLUSIONLIST)
        params:
            links = EXCLUSIONLIST_LINKS,
            outdir = 'exclusion_lists'
        log: f"logs/{TIMESTAMP}_download_exclusion_lists.log"
        threads: INIT_EXCLUSIONLIST_THREADS
        run:
            shell(f"echo '[{MODULE}] Downloading ENCODE blacklist exclusion regions...' | tee {{log}}")
            
            for i, link in enumerate(params.links, 1):
                filename = os.path.basename(link)
                new_filename = os.path.splitext(filename)[0]
                genome = new_filename.split('-')[0]
                
                shell(f"echo '[{MODULE}]   [{i}/{len(params.links)}] Downloading {genome}...' | tee -a {{log}}")
                shell(f"""
                    curl -L {link} --remote-name --silent
                    gunzip {filename}
                    sed -i 's/chr//g' {new_filename}
                    mv {new_filename} {{params.outdir}}/{new_filename}
                """)
            
            shell(f"echo '[{MODULE}] ✓ All exclusion lists downloaded ({len(params.links)} genomes)' | tee -a {{log}}")
            shell("touch {output.complete}")


# ================================================================================
# CHROMOSOME SIZES EXTRACTION
# ================================================================================

if INIT_CHRSIZESREF:
    rule get_chromosome_sizes:
        """
        Extracts chromosome sizes from reference genome FASTA file.
        Creates .sizes file required by bedtools genomecov and other tools
        for genome coverage calculations and bedGraph generation.
        """
        input:
            config['init']['chr_sizes_ref']
        output:
            INIT_CHRSIZESREF
        params:
            ref = os.path.basename(config['init']['chr_sizes_ref']),
            outdir = "chromosome_sizes"
        log: f"logs/{TIMESTAMP}_get_chromosome_sizes.log"
        threads: INIT_CHRSIZESREF_THREADS
        shell:
            f'''
            echo '[{MODULE}] Indexing reference genome...' | tee {{log}}
            samtools faidx {{input}} 2>> {{log}}
            
            echo '[{MODULE}] Extracting chromosome sizes...' | tee -a {{log}}
            awk '{{{{print $1"\\t"$2}}}}' {{input}}.fai > {{params.outdir}}/{{params.ref}}.sizes
            rm {{input}}.fai
            
            echo '[{MODULE}] ✓ Chromosome sizes extracted' | tee -a {{log}}
            '''


# ================================================================================
# KALLISTO INDEX BUILDING
# ================================================================================

if INIT_KALLISTO:
    rule kallisto_index:
        """
        Builds kallisto index for fast transcript quantification.
        Creates k-mer index from transcriptome FASTA for pseudoalignment-based
        quantification in RNA-seq workflows. Much faster than traditional alignment.
        """
        input:
            ref = config['init']['kallisto_ref']
        output:
            idx = INIT_KALLISTO
        log: f"logs/{TIMESTAMP}_kallisto_index.log"
        threads: INIT_KALLISTO_THREADS
        shell:
            f'''
            echo '[{MODULE}] Building kallisto k-mer index from transcriptome...' | tee {{log}}
            kallisto index --index {{output.idx}} {{input.ref}} 2>> {{log}}
            echo '[{MODULE}] ✓ Kallisto index complete' | tee -a {{log}}
            '''


# ================================================================================
# KRAKEN2/BRACKEN DATABASE SETUP
# ================================================================================

if INIT_KRAKEN2:
    rule kraken2_db:
        """
        Builds or downloads Kraken2 taxonomic classification database with Bracken support.
        
        Pre-built genomic databases (October 2025):
        All RefSeq databases include Bracken for 50, 75, 100, 150, 200, 250, and 300-mers.
        
        - viral: Refseq viral (0.6 GB)
        - minusb: No bacteria - archaea, viral, plasmid, human, UniVec_Core (11 GB)
        - standard: Full standard database - archaea, bacteria, viral, plasmid, human, UniVec_Core (94 GB)
        - standard_8/16: Memory-capped standard at 8/16 GB (7.5/14.9 GB)
        - pluspf: Standard + protozoa & fungi (100.6 GB)
        - pluspf_8/16: Memory-capped PlusPF at 8/16 GB (7.5/14.9 GB)
        - pluspfp: Standard + protozoa, fungi & plant (214.5 GB)
        - pluspfp_8/16: Memory-capped PlusPFP at 8/16 GB (7.5/14.9 GB)
        - pracken: Single reference per species - bacteria, archaea, protists, fungi (351.6 GB)
        - eupathdb48: Eukaryotic pathogens (11 GB)
        - core_nt: GenBank, RefSeq, TPA, PDB (316.2 GB)
        - gtdb_v226: GTDB v226 bacterial/archaeal genomes (644 GB)
        
        Pre-built 16S rRNA databases (March 2020):
        Include Bracken for 100, 150, and 200-mers.
        
        - 16s_greengenes_13.5: Greengenes 13.5 (73 MB)
        - 16s_rdp_11.5: RDP 11.5 (168 MB)
        - 16s_silva_132: Silva 132 (117 MB)
        - 16s_silva_138: Silva 138 (112 MB)
        
        Custom build types (from source):
        - greengenes/rdp/silva: Build 16S rRNA databases from scratch
        - nt_custom: NCBI nucleotide database (build from source, ~100GB+)
        - standard_custom: Build standard database from source
        
        Note: 
        - RefSeq databases use default k-mer (35), minimizer length, and spacing
        - Capped databases (-8, -16) use --max-db-size for reduced sensitivity/accuracy
        - All indexes include inspect.txt, ktaxonomy.tsv, library_report.tsv, and taxonomy files
        - PrackenDB is optimized for k-mer counting (single genome per species)
        - Custom builds generate Bracken for multiple read lengths
        """
        output:
            db = INIT_KRAKEN2
        params:
            db_type = config['init']['kraken2_db_type'],
            read_len = config.get('init', {}).get('bracken_read_len', 100),
            outdir = 'kraken2_db/' + config['init']['kraken2_db_type']
        log: f"logs/{TIMESTAMP}_kraken2_db.log"
        threads: INIT_KRAKEN2_THREADS
        run:
            shell(f"echo '[{MODULE}] Setting up Kraken2/Bracken database: {{params.db_type}}...' | tee {{log}}")

            # Check if this is a pre-built database
            if params['db_type'] in KRAKEN2_PREBUILT_DBS:
                url, archive_size, index_size, description = KRAKEN2_PREBUILT_DBS[params['db_type']]
                filename = os.path.basename(url)
                
                # Special handling for 16S databases
                if params['db_type'].startswith('16s_'):
                    shell(f"""
                        echo '[{MODULE}] Downloading pre-built {description}...' | tee -a {{log}}
                        echo '[{MODULE}] Size: {archive_size} GB (includes Bracken 100/150/200mers)' | tee -a {{log}}
                        curl -L {url} --remote-name --silent --show-error 2>> {{log}}
                        
                        echo '[{MODULE}] Extracting database...' | tee -a {{log}}
                        tar -zxf {filename} -C {{params.outdir}} 2>> {{log}}
                        rm {filename}
                        
                        echo '[{MODULE}] ✓ 16S rRNA database ready with Bracken support' | tee -a {{log}}
                        touch {{output.db}}
                    """)
                
                # Special handling for PrackenDB
                elif params['db_type'] == 'pracken':
                    shell(f"""
                        echo '[{MODULE}] Downloading PrackenDB (single genome per species)...' | tee -a {{log}}
                        echo '[{MODULE}] Archive size: {archive_size} GB' | tee -a {{log}}
                        echo '[{MODULE}] Optimized for k-mer counting methods' | tee -a {{log}}
                        curl -L {url} --remote-name --silent --show-error 2>> {{log}}
                        
                        echo '[{MODULE}] Extracting database...' | tee -a {{log}}
                        tar -zxf {filename} -C {{params.outdir}} 2>> {{log}}
                        rm {filename}
                        
                        echo '[{MODULE}] ✓ PrackenDB downloaded and extracted' | tee -a {{log}}
                        touch {{output.db}}
                    """)
                
                # Standard pre-built RefSeq databases (with Bracken 50/75/100/150/200/250/300mers)
                else:
                    shell(f"""
                        echo '[{MODULE}] Downloading pre-built {description}...' | tee -a {{log}}
                        echo '[{MODULE}] Archive size: {archive_size} GB, Index size: {index_size} GB' | tee -a {{log}}
                        echo '[{MODULE}] Includes Bracken for 50, 75, 100, 150, 200, 250, 300-mers' | tee -a {{log}}
                        curl -L {url} --remote-name --silent --show-error 2>> {{log}}
                        
                        echo '[{MODULE}] Extracting database...' | tee -a {{log}}
                        tar -zxf {filename} -C {{params.outdir}} 2>> {{log}}
                        rm {filename}
                        
                        echo '[{MODULE}] ✓ Pre-built database downloaded and extracted' | tee -a {{log}}
                        echo '[{MODULE}] Database includes: inspect.txt, ktaxonomy.tsv, library_report.tsv' | tee -a {{log}}
                        touch {{output.db}}
                    """)
            
            # Custom build: Standard RefSeq database
            elif params['db_type'] == 'standard_custom':
                shell(f"""
                    echo '[{MODULE}] Building standard database from source (bacteria/archaea/viral/human/plasmid)...' | tee -a {{log}}
                    kraken2-build --threads {{threads}} --use-ftp --standard --db {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] Building Bracken for {{params.read_len}}mers...' | tee -a {{log}}
                    bracken-build -t {{threads}} -k 35 -l {{params.read_len}} -d {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] ✓ Standard custom database complete (Bracken {{params.read_len}}mers)' | tee -a {{log}}
                """)

            # Custom build: Specialized 16S rRNA databases
            elif params['db_type'] in ['greengenes','rdp','silva']:
                shell(f"""
                    echo '[{MODULE}] Building {{params.db_type}} 16S rRNA database from source...' | tee -a {{log}}
                    kraken2-build --threads {{threads}} --use-ftp --special {{params.db_type}} --db {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] Building Bracken for {{params.read_len}}mers...' | tee -a {{log}}
                    bracken-build -t {{threads}} -k 35 -l {{params.read_len}} -d {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] ✓ 16S custom database complete (Bracken {{params.read_len}}mers)' | tee -a {{log}}
                """)

            # Custom build: NCBI nucleotide database (very large)
            elif params['db_type'] == 'nt_custom':
                shell(f"""
                    echo '[{MODULE}] Downloading NCBI taxonomy...' | tee -a {{log}}
                    kraken2-build --use-ftp --download-taxonomy --db {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] Downloading NCBI nucleotide library (this may take hours)...' | tee -a {{log}}
                    kraken2-build --use-ftp --download-library nt --db {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] Building Kraken2 database...' | tee -a {{log}}
                    kraken2-build --threads {{threads}} --build --fast-build --db {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] Building Bracken for {{params.read_len}}mers...' | tee -a {{log}}
                    bracken-build -t {{threads}} -k 35 -l {{params.read_len}} -d {{params.outdir}} 2>> {{log}}
                    
                    echo '[{MODULE}] ✓ NCBI nucleotide custom database complete (Bracken {{params.read_len}}mers)' | tee -a {{log}}
                """)

            else:
                shell(f"""
                    echo '[{MODULE}] ✗ ERROR: Unknown database type {{params.db_type}}' | tee -a {{log}}
                    echo '[{MODULE}] Available pre-built genomic databases:' | tee -a {{log}}
                    echo '[{MODULE}]   viral, minusb, standard, standard_8, standard_16,' | tee -a {{log}}
                    echo '[{MODULE}]   pluspf, pluspf_8, pluspf_16, pluspfp, pluspfp_8, pluspfp_16,' | tee -a {{log}}
                    echo '[{MODULE}]   pracken, eupathdb48, core_nt, gtdb_v226' | tee -a {{log}}
                    echo '[{MODULE}]   (RefSeq databases include Bracken 50/75/100/150/200/250/300mers)' | tee -a {{log}}
                    echo '[{MODULE}] Available pre-built 16S databases:' | tee -a {{log}}
                    echo '[{MODULE}]   16s_greengenes_13.5, 16s_rdp_11.5, 16s_silva_132, 16s_silva_138' | tee -a {{log}}
                    echo '[{MODULE}]   (16S databases include Bracken 100/150/200mers)' | tee -a {{log}}
                    echo '[{MODULE}] Custom build types:' | tee -a {{log}}
                    echo '[{MODULE}]   standard_custom, greengenes, rdp, silva, nt_custom' | tee -a {{log}}
                    exit 1
                """)

            shell(f"echo '[{MODULE}] ✓ Kraken2/Bracken database setup complete' | tee -a {{log}}")


# ================================================================================
# WORKFLOW COMPLETION
# ================================================================================

onsuccess:
    # Build list of completed components
    completed = []
    output_dirs = []
    
    if INIT_HISAT2:
        gtf_status = "with GTF annotation" if config["init"].get("hisat2_gtf") else "no GTF"
        completed.append(f"HISAT2 index ({gtf_status})")
        output_dirs.append(f"hisat2_index/ (alignment index)")
    
    if INIT_KALLISTO:
        completed.append("Kallisto index")
        output_dirs.append("kallisto_index/ (transcript quantification)")
    
    if INIT_KRAKEN2:
        db_type = config['init']['kraken2_db_type']
        if db_type in KRAKEN2_PREBUILT_DBS:
            completed.append(f"Kraken2/Bracken: {db_type} (pre-built)")
        else:
            completed.append(f"Kraken2/Bracken: {db_type} (custom)")
        output_dirs.append(f"kraken2_db/{db_type}/ (taxonomic classification)")
    
    if INIT_RIBORNA:
        completed.append("rRNA references (SortMeRNA)")
        output_dirs.append("rRNA_ref/ (contamination removal)")
    
    if INIT_EXCLUSIONLIST:
        completed.append(f"Exclusion lists ({len(EXCLUSIONLIST_LINKS)} genomes)")
        output_dirs.append("exclusion_lists/ (blacklisted regions)")
    
    if INIT_CHRSIZESREF:
        completed.append("Chromosome sizes")
        output_dirs.append("chromosome_sizes/ (genome coverage)")
    
    completed_str = "\n".join([f"[{MODULE}]   • {c}" for c in completed])
    dirs_str = "\n".join([f"[{MODULE}]   - {d}" for d in output_dirs])
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✓ Database initialization complete!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Initialized {len(completed)} component(s):'
        echo '{completed_str}'
        echo ''
        echo '[{MODULE}] Output directories:'
        echo '{dirs_str}'
        echo ''
        echo '[{MODULE}] All reference files ready for analysis workflows'
        echo ''
    """)


onerror:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✗ Database initialization FAILED!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Common issues:'
        echo '[{MODULE}]   - Insufficient disk space for databases'
        echo '[{MODULE}]   - Network errors during downloads'
        echo '[{MODULE}]   - Missing tools (hisat2-build, kallisto, kraken2, etc.)'
        echo '[{MODULE}]   - Invalid reference paths or file formats'
        echo '[{MODULE}]   - Memory issues during index building'
        echo '[{MODULE}]   - Timeout on large database downloads'
        echo ''
        echo '[{MODULE}] Troubleshooting:'
        echo '[{MODULE}]   1. Check logs/ directory for detailed errors'
        echo '[{MODULE}]   2. Verify sufficient disk space for target databases'
        echo '[{MODULE}]   3. Confirm all required tools are in PATH'
        echo '[{MODULE}]   4. Test network connectivity for downloads'
        echo '[{MODULE}]   5. Validate input file paths in configuration'
        echo ''
        echo '[{MODULE}] Failed at timestamp: {TIMESTAMP}'
        echo ''
    """)
