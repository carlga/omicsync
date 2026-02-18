import os
import sys
import re
import glob
import subprocess

from .logger import logger


#
# init module
#
def init_config_set(config, args):
    """ Sets parameters for init module. """

    if config['init']['hisat2_ref'] is None:
        if args.hisat2_ref:
            config['init']['hisat2_ref'] = os.path.abspath(args.hisat2_ref)
    else:
        logger.warning(
            f"HISAT2 reference file '{config['init']['hisat2_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['hisat2_index_base'] is None:
        config['init']['hisat2_index_base'] = args.hisat2_index_base
    else:
        logger.warning(
            f"HISAT2 index basename '{config['init']['hisat2_index_base']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['hisat2_gtf'] is None:
        if args.hisat2_gtf:
            config['init']['hisat2_gtf'] = os.path.abspath(args.hisat2_gtf)
    else:
        logger.warning(
            f"HISAT2 annotation file '{config['init']['hisat2_gtf']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['rrna_ref'] is None:
        config['init']['rrna_ref'] = args.rrna_ref
    else:
        logger.warning(
            f"Download rRNA reference '{config['init']['rrna_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['exclusion_lists'] is None:
        config['init']['exclusion_lists'] = args.exclusion_lists
    else:
        logger.warning(
            f"Download exclusion lists '{config['init']['exclusion_lists']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['chr_sizes_ref'] is None:
        if args.chr_sizes_ref:
            config['init']['chr_sizes_ref'] = os.path.abspath(args.chr_sizes_ref)
    else:
        logger.warning(
            f"Chromosome sizes reference file '{config['init']['chr_sizes_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['kallisto_ref'] is None: 
        if args.kallisto_ref:
            config['init']['kallisto_ref'] = os.path.abspath(args.kallisto_ref)
    else:
        logger.warning(
            f"Kallisto reference file '{config['init']['kallisto_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['kallisto_index'] is None:
        config['init']['kallisto_index'] = args.kallisto_index
    else:
        logger.warning(
            f"Kallisto index filename '{config['init']['kallisto_index']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['kraken2_db_type'] is None:
        config['init']['kraken2_db_type'] = args.kraken2_db_type
    else:
        logger.warning(
            f"Kraken2 database '{config['init']['kraken2_db_type']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['bracken_read_len'] is None:
        if config['init']['kraken2_db_type']:
            config['init']['bracken_read_len'] = args.bracken_read_len
    else:
        logger.warning(
            f"Bracken read length '{config['init']['bracken_read_len']}' "
            "is set in configuration! Using..."
        )
    
    if config['init']['outdir'] is None:
        config['init']['outdir'] = os.path.abspath(args.outdir)
    else:
        logger.warning(
            f"Outdir '{config['init']['outdir']}' is set in configuration! Using..."
        )

    if config['init']['jobs'] is None:
        config['init']['jobs'] = args.jobs
    else:
        logger.warning(
            f"Job/core number '{config['init']['jobs']}' is set in configuration! Using..."
        )
    
    if config['init']['profile'] is None:
        if args.profile:
            config['init']['profile'] = os.path.abspath(args.profile)
    else:
        logger.warning(
            f"Profile '{config['init']['profile']}' is set in configuration! Using..."
        )
    
    return config


def init_config_check(config):
    """ Evaluates input parameters for init module. """

    # HISAT2 index build
    if config['init']['hisat2_ref']:
        if not os.path.exists(config['init']['hisat2_ref']):
            logger.error(f"HISAT2 reference '{config['init']['hisat2_ref']}' not found.")
            sys.exit(1)
    
    if config['init']['hisat2_ref'] and config['init']['hisat2_index_base'] is None:
        logger.error("Missing basename for HISAT2 index.")
        sys.exit(1)

    if config['init']['hisat2_gtf']:
        if not os.path.exists(config['init']['hisat2_gtf']):
            logger.error(f"Annotation file '{config['init']['hisat2_gtf']}' not found.")
            sys.exit(1)
    
    # Reference genome for chromosome sizes
    if config['init']['chr_sizes_ref']:
        if not os.path.exists(config['init']['chr_sizes_ref']):
            logger.error(f"Chromosome sizes reference '{config['init']['chr_sizes_ref']}' not found.")
            sys.exit(1)

    # kallisto index build
    if config['init']['kallisto_ref']:
        if not os.path.exists(config['init']['kallisto_ref']):
            logger.error(f"kallisto reference '{config['init']['kallisto_ref']}' not found.")
            sys.exit(1)
    
    if config['init']['kallisto_ref'] and config['init']['kallisto_index'] is None:
        logger.error("Missing filename for kallisto index.")
        sys.exit(1)

    # kraken2/bracken DB build
    avail_k2db = [
        'standard', 'greengenes', 'rdp', 'silva', 'nt', 'std16_prebuilt', 
        'plusPFP16_prebuilt'
    ]
    if config['init']['kraken2_db_type']:
        if config['init']['kraken2_db_type'] not in avail_k2db:
            logger.error(
                "kraken2 DB is not valid! Choose from: 'standard', 'greengenes', "
                "'rdp', 'silva', 'nt', 'std16_prebuilt' and 'plusPFP16_prebuilt'.")
            sys.exit(1)
    
    # Profile
    if config['init']['profile']:
        if not os.path.exists(config['init']['profile']):
            logger.error(f"Snakemake profile '{config['init']['profile']}' not found.")
            sys.exit(1)

    
#
# getref module
#
def getref_config_set(config, args):
    """ Sets parameters for getref module. """

    if config['getref']['release'] is None:
        # Set latest Ensembl release if none specified
        cmd = "gget ref homo_sapiens -w dna"
        ensembl_latest = subprocess.check_output(
            cmd, 
            shell=True, 
            text=True, 
            stderr=subprocess.DEVNULL
        )
        ensembl_latest = int(re.search(r'"ensembl_release":\s+(\d+),', ensembl_latest).group(1))
        config['getref']['release'] = args.release if args.release is not None else ensembl_latest
    else:
        logger.warning(
            f"Release '{config['getref']['release']}' is set in configuration! Using..."
        )
    
    # List available species
    if args.list_sp:
        subprocess.check_call(
            f"gget ref --list_species --release {config['getref']['release']}", 
            shell=True
        )
        sys.exit(0)
    
    if config['getref']['species'] is None:
        config['getref']['species'] = args.species
    else:
        logger.warning(
            f"Species '{config['getref']['species']}' is set in configuration! Using..."
        )
    
    if config['getref']['which'] is None:
        config['getref']['which'] = args.which
    else:
        logger.warning(
            f"Reference files '{config['getref']['which']}' are set in configuration! Using..."
        )
    
    if config['getref']['outdir'] is None:
        config['getref']['outdir'] = os.path.abspath(args.outdir)
    else:
        logger.warning(
            f"Outdir '{config['getref']['outdir']}' is set in configuration! Using..."
        )
    
    if config['getref']['jobs'] is None:
        config['getref']['jobs'] = args.jobs
    else:
        logger.warning(
            f"Job/core number '{config['getref']['jobs']}' is set in configuration! Using..."
        )
    
    if config['getref']['profile'] is None:
        if args.profile:
            config['getref']['profile'] = os.path.abspath(args.profile)
    else:
        logger.warning(
            f"Profile '{config['getref']['profile']}' is set in configuration! Using..."
        )
    
    return config


def getref_config_check(config):
    """ Evaluates input parameters for getref module. """

    # Ensembl release
    cmd = "gget ref homo_sapiens -w dna"
    ensembl_latest = subprocess.check_output(
        cmd, 
        shell=True, 
        text=True, 
        stderr=subprocess.DEVNULL
    )
    ensembl_latest = int(re.search(r'"ensembl_release":\s+(\d+),', ensembl_latest).group(1))
    if config['getref']['release'] > ensembl_latest:
        logger.error(f"Ensembl release '{config['getref']['release']}' is not valid.")
        sys.exit(1)
    
    # Species
    if config['getref']['species'] is None:
        logger.error("Species is not specified.")
        sys.exit(1)

    # Profile
    if config['getref']['profile']:
        if not os.path.exists(config['getref']['profile']):
            logger.error(f"Snakemake profile '{config['getref']['profile']}' not found.")
            sys.exit(1)


#
# getfq module
#
def getfq_config_set(config, args):
    """ Sets parameters for getfq module. """

    if 'getfq' not in config:
        config['getfq'] = {}

    # Mutually exclusive ID sources
    has_cmdline = hasattr(args, 'ids') and args.ids
    has_file = hasattr(args, 'ids_file') and args.ids_file
    
    if has_cmdline and has_file:
        logger.error("Cannot use both command-line IDs and --ids-file simultaneously")
        sys.exit(1)
    
    # Load IDs
    all_ids = []
    if has_file:
        if not os.path.exists(args.ids_file):
            logger.error(f"IDs file not found: {args.ids_file}")
            sys.exit(1)
        with open(args.ids_file, 'r') as f:
            all_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        logger.warning(f"Loaded {len(all_ids)} IDs from file")
    elif has_cmdline:
        all_ids = args.ids
    else:
        logger.error("No accession IDs provided (use command-line or --ids-file)")
        sys.exit(1)
    
    # Deduplicate
    unique_ids = list(dict.fromkeys(all_ids))
    if len(unique_ids) < len(all_ids):
        logger.warning(f"Removed {len(all_ids) - len(unique_ids)} duplicates")
    
    config['getfq']['ids'] = ' '.join(unique_ids)
    logger.warning(f"Total IDs: {len(unique_ids)}")

    # Other settings
    if hasattr(args, 'ftp'): config['getfq']['ftp'] = args.ftp
    if hasattr(args, 'batch_size'): 
        config['getfq']['batch_size'] = args.batch_size
        if args.batch_size > 1: logger.warning(f"Batch size: {args.batch_size}")
    if hasattr(args, 'outdir') and args.outdir: config['getfq']['outdir'] = args.outdir
    if hasattr(args, 'jobs') and args.jobs: config['getfq']['jobs'] = args.jobs
    if hasattr(args, 'profile'): config['getfq']['profile'] = args.profile
    if hasattr(args, 'dryrun'): config['getfq']['dryrun'] = args.dryrun

    return config


def getfq_config_check(config):
    """ Evaluates input parameters for getfq module. """

    if 'getfq' not in config:
        logger.error("getfq configuration not found")
        sys.exit(1)

    if 'ids' not in config['getfq'] or not config['getfq']['ids']:
        logger.error("No accession IDs provided")
        sys.exit(1)

    if 'outdir' in config['getfq']:
        outdir = config['getfq']['outdir']
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logger.warning(f"Created output directory: {outdir}")
    
    batch_size = config['getfq'].get('batch_size', 1)
    if batch_size < 1:
        logger.error(f"Invalid batch size: {batch_size}")
        sys.exit(1)
    
    # Profile
    if config['getfq']['profile']:
        if not os.path.exists(config['getfq']['profile']):
            logger.error(f"Snakemake profile '{config['getfq']['profile']}' not found.")
            sys.exit(1)


#
# runqc module
#
def runqc_config_set(config, args):
    """ Sets parameters for runqc module. """

    if config['runqc']['indir'] is None:
        if args.indir:
            config['runqc']['indir'] = os.path.abspath(args.indir)
    else:
        logger.warning(
            f"Input directory '{config['runqc']['indir']}' is set in configuration! Using..."
        )
    
    if config['runqc']['number'] is None:
        # Set to integer when number > 1
        number = int(args.number) if args.number > 1 else args.number
        config['runqc']['number'] = number
    else:
        logger.warning(
            f"Sampling number '{config['runqc']['number']}' is set in configuration! Using..."
        )
    
    if config['runqc']['seed'] is None:
        config['runqc']['seed'] = args.seed
    else:
        logger.warning(
            f"Random seed '{config['runqc']['seed']}' is set in configuration! Using..."
        )
    
    if config['runqc']['fastp_opt'] is None:
        config['runqc']['fastp_opt'] = args.fastp_opt
    else:
        logger.warning(
            f"fastp options '{config['runqc']['fastp_opt']}' are set in configuration! Using..."
        )
    
    if config['runqc']['hisat2_ref'] is None:
        if args.hisat2_ref:
            config['runqc']['hisat2_ref'] = os.path.abspath(args.hisat2_ref)
    else:
        logger.warning(
            f"HISAT2 index '{config['runqc']['hisat2_ref']}' is set in configuration! Using..."
        )
    
    if config['runqc']['gtf'] is None:
        if args.gtf:
            config['runqc']['gtf'] = os.path.abspath(args.gtf)
    else:
        logger.warning(
            f"Annotation file '{config['runqc']['gtf']}' is set in configuration! Using..."
        )
    
    if config['runqc']['kraken2_db'] is None:
        if args.kraken2_db:
            config['runqc']['kraken2_db'] = os.path.abspath(args.kraken2_db) 
    else:
        logger.warning(
            f"Kraken DB '{config['runqc']['kraken2_db']}' is set in configuration! Using..."
        )
    
    if config['runqc']['outdir'] is None:
        config['runqc']['outdir'] = os.path.abspath(args.outdir)
    else:
        logger.warning(
            f"Outdir '{config['runqc']['outdir']}' is set in configuration! Using..."
        )
    
    if config['runqc']['jobs'] is None:
        config['runqc']['jobs'] = args.jobs
    else:
        logger.warning(
            f"Job/core number '{config['runqc']['jobs']}' is set in configuration! Using..."
        )
    
    if config['runqc']['profile'] is None:
        if args.profile:
            config['runqc']['profile'] = os.path.abspath(args.profile)
    else:
        logger.warning(
            f"Profile '{config['runqc']['profile']}' is set in configuration! Using..."
        )
    
    return config


def runqc_config_check(config):
    """ Evaluates input parameters for runqc module. """

    # Input directory
    if config['runqc']['indir'] is None:
        logger.error("Missing fastq file input directory.")
        sys.exit(1)
    
    if config['runqc']['indir']:
        if not os.path.exists(config['runqc']['indir']):
            logger.error(f"Input directory '{config['runqc']['indir']}' not found.")
            sys.exit(1)
        elif not glob.glob(os.path.join(config['runqc']['indir'], '*.f*q.gz')):
            logger.error(f"No fastq files found in '{config['runqc']['indir']}'.")
            sys.exit(1)

    # fastp options
    blocked_opts = [
        '-i', '--in1', '-o', '--out1', '-I', '--in2', '-O', '--out2', '-m', '--merge',
        '--merged_out', '--include_unmerged', '--stdin', '--stdout', '--interleaved_in', 
        '-j', '--json', '-h', '--html', '-s', '--split', '-S', '--split_by_lines', 
        '-d', '--split_prefix_digits', '-?', '--help', '-v', '--version'
    ]
    if config['runqc']['fastp_opt']:
        fastp_opts = re.split('\s+', config['runqc']['fastp_opt'])
        for opt in fastp_opts:
            if opt in blocked_opts:
                logger.error(f"fastp option '{opt}' is not allowed.")
                sys.exit(1)
    
    # HISAT2 index
    if config['runqc']['hisat2_ref']:
        if not os.path.exists(config['runqc']['hisat2_ref'] + '.1.ht2'):
            logger.error(f"Missing HISAT2 index files at '{config['runqc']['hisat2_ref']}'.")
            sys.exit(1)
    
    # Annotation
    if config['runqc']['gtf']:
        if not os.path.exists(config['runqc']['gtf']):
            logger.error(f"Annotation file '{config['runqc']['gtf']}' not found.")
            sys.exit(1)
    
    # kraken2/bracken DB
    if config['runqc']['kraken2_db']:
        # kraken2 
        if not os.path.exists(config['runqc']['kraken2_db']):
            logger.error(f"kraken2/bracken DB '{config['runqc']['kraken2_db']}' not found.")
            sys.exit(1)
        elif not os.path.exists(os.path.join(config['runqc']['kraken2_db'], 'hash.k2d')):
            logger.error(f"Missing 'hash.k2d' in kraken2 DB.")
            sys.exit(1)
        elif not os.path.exists(os.path.join(config['runqc']['kraken2_db'], 'opts.k2d')):
            logger.error(f"Missing 'opts.k2d' in kraken2 DB.")
            sys.exit(1)
        elif not os.path.exists(os.path.join(config['runqc']['kraken2_db'], 'seqid2taxid.map')):
            logger.error(f"Missing 'seqid2taxid.map' in kraken2 DB.")
            sys.exit(1)
        elif not os.path.exists(os.path.join(config['runqc']['kraken2_db'], 'taxo.k2d')):
            logger.error(f"Missing 'taxo.k2d' in kraken2 DB.")
            sys.exit(1)
        
        # bracken
        kmer_distribs = glob.glob(os.path.join(config['runqc']['kraken2_db'], '*.kmer_distrib'))
        if not kmer_distribs:
            logger.error(f"No '.kmer_distrib' files found in '{config['runqc']['kraken2_db']}'.")
            sys.exit(1)
        if len(kmer_distribs) > 1:
            kmer_distr100 = os.path.join(config['runqc']['kraken2_db'], 'database100mers.kmer_distrib')
            if not os.path.exists(kmer_distr100):
                logger.error(
                    f"Multiple '.kmer_distrib' files found in '{config['runqc']['kraken2_db']}' "
                    "but missing 'database100mers.kmer_distrib'. To use a different kmer length "
                    "leave only one file in DB directory."
                )
                sys.exit(1)
            else:
                logger.warning(
                    f"Multiple '.kmer_distrib' files found in '{config['runqc']['kraken2_db']}'. "
                    "Using 'database100mers.kmer_distrib'..."
                )

    # Profile
    if config['runqc']['profile']:
        if not os.path.exists(config['runqc']['profile']):
            logger.error(f"Snakemake profile '{config['runqc']['profile']}' not found.")
            sys.exit(1)


#
# run_rnaseq module
#
def run_rnaseq_config_set(config, args):
    """ Sets parameters for run_rnaseq module. """

    if config['run_rnaseq']['indir'] is None:
        if args.indir:
            config['run_rnaseq']['indir'] = os.path.abspath(args.indir)
    else:
        logger.warning(
            f"Input directory '{config['run_rnaseq']['indir']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_rnaseq']['fastp_opt'] is None:
        config['run_rnaseq']['fastp_opt'] = args.fastp_opt
    else:
        logger.warning(
            f"fastp options '{config['run_rnaseq']['fastp_opt']}' "
            "are set in configuration! Using..."
        )

    if config['run_rnaseq']['removal_ref'] is None:
        if args.removal_ref:
            config['run_rnaseq']['removal_ref'] = os.path.abspath(args.removal_ref)
    else:
        logger.warning(
            f"Read removal reference '{config['run_rnaseq']['removal_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_rnaseq']['hisat2_ref'] is None:
        if args.hisat2_ref:
            config['run_rnaseq']['hisat2_ref'] = os.path.abspath(args.hisat2_ref)
    else:
        logger.warning(
            f"HISAT2 index '{config['run_rnaseq']['hisat2_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_rnaseq']['hisat2_opt'] is None:
        config['run_rnaseq']['hisat2_opt'] = args.hisat2_opt
    else:
        logger.warning(
            f"HISAT2 options '{config['run_rnaseq']['hisat2_opt']}' "
            "are set in configuration! Using..."
        )
    
    if config['run_rnaseq']['kallisto_ref'] is None:
        if args.kallisto_ref:
            config['run_rnaseq']['kallisto_ref'] = os.path.abspath(args.kallisto_ref)
    else:
        logger.warning(
            f"kallisto index '{config['run_rnaseq']['kallisto_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_rnaseq']['kallisto_opt'] is None:
        config['run_rnaseq']['kallisto_opt'] = args.kallisto_opt
    else:
        logger.warning(
            f"kallisto options '{config['run_rnaseq']['kallisto_opt']}' "
            "are set in configuration! Using..."
        )
    
    if config['run_rnaseq']['gtf'] is None:
        if args.gtf:
            config['run_rnaseq']['gtf'] = os.path.abspath(args.gtf)
    else:
        logger.warning(
            f"Annotation file '{config['run_rnaseq']['gtf']}' "
            "is set in configuration! Using..."
        )

    if config['run_rnaseq']['exon_quant'] is None:
        if args.exon_quant:
            config['run_rnaseq']['exon_quant'] = args.exon_quant
    else:
        logger.warning(
            f"Exon-level quantification '{config['run_rnaseq']['exon_quant']}' "
            "is set in configuration! Using..."
        )

    if config['run_rnaseq']['genome'] is None:
        if args.genome:
            config['run_rnaseq']['genome'] = os.path.abspath(args.genome)
    else:
        logger.warning(
            f"Reference FASTA file '{config['run_rnaseq']['genome']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_rnaseq']['outdir'] is None:
        config['run_rnaseq']['outdir'] = os.path.abspath(args.outdir)
    else:
        logger.warning(
            f"Outdir '{config['run_rnaseq']['outdir']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_rnaseq']['jobs'] is None:
        config['run_rnaseq']['jobs'] = args.jobs
    else:
        logger.warning(
            f"Job/core number '{config['run_rnaseq']['jobs']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_rnaseq']['profile'] is None:
        if args.profile:
            config['run_rnaseq']['profile'] = os.path.abspath(args.profile)
    else:
        logger.warning(
            f"Profile '{config['run_rnaseq']['profile']}' "
            "is set in configuration! Using..."
        )
    
    return config


def run_rnaseq_config_check(config):
    """ Evaluates input parameters for run_rnaseq module. """

    # Input directory
    if config['run_rnaseq']['indir'] is None:
        logger.error("Missing fastq file input directory.")
        sys.exit(1)
    else:
        if not os.path.exists(config['run_rnaseq']['indir']):
            logger.error(f"Input directory '{config['run_rnaseq']['indir']}' not found.")
            sys.exit(1)
        elif not glob.glob(os.path.join(config['run_rnaseq']['indir'], '*.f*q.gz')):
            logger.error(f"No fastq files found in '{config['run_rnaseq']['indir']}'.")
            sys.exit(1)

    # fastp options
    blocked_opts = [
        '-i', '--in1', '-o', '--out1', '-I', '--in2', '-O', '--out2', '-m', '--merge',
        '--merged_out', '--include_unmerged', '--stdin', '--stdout', '--interleaved_in', 
        '-j', '--json', '-h', '--html', '-s', '--split', '-S', '--split_by_lines', 
        '-d', '--split_prefix_digits', '-?', '--help', '-v', '--version'
    ]
    if config['run_rnaseq']['fastp_opt']:
        fastp_opts = re.split('\s+', config['run_rnaseq']['fastp_opt'])
        for opt in fastp_opts:
            if opt in blocked_opts:
                logger.error(f"fastp option '{opt}' is not allowed.")
                sys.exit(1)
    
    # Read removal reference
    if config['run_rnaseq']['removal_ref']:
        if not os.path.exists(config['run_rnaseq']['removal_ref']):
            logger.error(f"Read removal reference file '{config['run_rnaseq']['removal_ref']}' not found.")
            sys.exit(1)

    # HISAT2 index and options
    if config['run_rnaseq']['hisat2_ref'] is None:
        logger.error("Missing HISAT2 index.")
        sys.exit(1)
    else:
        if not os.path.exists(config['run_rnaseq']['hisat2_ref'] + '.1.ht2'):
            logger.error(f"Missing HISAT2 index files at '{config['run_rnaseq']['hisat2_ref']}'.")
            sys.exit(1)
    
    blocked_opts = [
        '-x', '-1', '-2', '-U', '-S', '--qseq', '-f', '-r', '-c', '--nofw', '--norc', 
        '--summary-file', '--version', '-h', '--help'
    ]
    if config['run_rnaseq']['hisat2_opt']:
        hisat2_opts = re.split('\s+', config['run_rnaseq']['hisat2_opt'])
        for opt in hisat2_opts:
            if opt in blocked_opts:
                logger.error(f"HISAT2 option '{opt}' is not allowed.")
                sys.exit(1)
    
    # kallisto index and options
    if config['run_rnaseq']['kallisto_ref']:
        if not os.path.exists(config['run_rnaseq']['kallisto_ref']):
            logger.error(f"kallisto index file '{config['run_rnaseq']['kallisto_ref']}' not found.")
            sys.exit(1)
    
    blocked_opts = [
        '-i', '--index', '-o', '--output-dir', '--plaintext', '--single'
    ]
    if config['run_rnaseq']['kallisto_opt']:
        kallisto_opts = re.split('\s+', config['run_rnaseq']['kallisto_opt'])
        for opt in kallisto_opts:
            if opt in blocked_opts:
                logger.error(f"kallisto option '{opt}' is not allowed.")
                sys.exit(1)

    # Annotation
    if config['run_rnaseq']['gtf'] is None:
        logger.error("Missing gene annotation file.")
        sys.exit(1)
    else:
        if not os.path.exists(config['run_rnaseq']['gtf']):
            logger.error(f"Annotation file '{config['run_rnaseq']['gtf']}' not found.")
            sys.exit(1)
    
    # Exon-level quantification
    if config['run_rnaseq']['exon_quant']:
        if config['run_rnaseq']['genome']:
            if not os.path.exists(config['run_rnaseq']['genome']):
                logger.error(f"Reference FASTA file '{config['run_rnaseq']['genome']}' not found.")
                sys.exit(1)
    else:
        if config['run_rnaseq']['genome']:
            logger.error(f"Reference FASTA file '{config['run_rnaseq']['genome']}' should be provided with -e/--exon_quent flag for exon-level quantification.")
            sys.exit(1)
    
    # Profile
    if config['run_rnaseq']['profile']:
        if not os.path.exists(config['run_rnaseq']['profile']):
            logger.error(f"Snakemake profile '{config['run_rnaseq']['profile']}' not found.")
            sys.exit(1)


#
# run_atacseq module
#
def run_atacseq_config_set(config, args):
    """ Sets parameters for run_atacseq module. """

    if config['run_atacseq']['indir'] is None:
        if args.indir:
            config['run_atacseq']['indir'] = os.path.abspath(args.indir)
    else:
        logger.warning(
            f"Input directory '{config['run_atacseq']['indir']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['fastp_opt'] is None:
        config['run_atacseq']['fastp_opt'] = args.fastp_opt
    else:
        logger.warning(
            f"fastp options '{config['run_atacseq']['fastp_opt']}' "
            "are set in configuration! Using..."
        )
    
    if config['run_atacseq']['hisat2_ref'] is None:
        if args.hisat2_ref:
            config['run_atacseq']['hisat2_ref'] = os.path.abspath(args.hisat2_ref)
    else:
        logger.warning(
            f"HISAT2 index '{config['run_atacseq']['hisat2_ref']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['hisat2_opt'] is None:
        config['run_atacseq']['hisat2_opt'] = args.hisat2_opt
    else:
        logger.warning(
            f"HISAT2 options '{config['run_atacseq']['hisat2_opt']}' "
            "are set in configuration! Using..."
        )
    
    if config['run_atacseq']['exclusion_list'] is None:
        if args.exclusion_list:
            config['run_atacseq']['exclusion_list'] = os.path.abspath(args.exclusion_list)
    else:
        logger.warning(
            f"Exclusion list '{config['run_atacseq']['exclusion_list']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['chr_sizes'] is None:
        if args.chr_sizes:
            config['run_atacseq']['chr_sizes'] = os.path.abspath(args.chr_sizes)
    else:
        logger.warning(
            f"Chromosome size file '{config['run_atacseq']['chr_sizes']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['tn5_shift'] is None:
        config['run_atacseq']['tn5_shift'] = args.tn5_shift
    else:
        logger.warning(
            f"Tn5 shift '{config['run_atacseq']['tn5_shift']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['genome_size'] is None:
        config['run_atacseq']['genome_size'] = args.genome_size
    else:
        logger.warning(
            f"Effective genome size '{config['run_atacseq']['genome_size']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['gtf'] is None:
        if args.gtf:
            config['run_atacseq']['gtf'] = os.path.abspath(args.gtf)
    else:
        logger.warning(
            f"Annotation file '{config['run_atacseq']['gtf']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['outdir'] is None:
        config['run_atacseq']['outdir'] = os.path.abspath(args.outdir)
    else:
        logger.warning(
            f"Outdir '{config['run_atacseq']['outdir']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['jobs'] is None:
        config['run_atacseq']['jobs'] = args.jobs
    else:
        logger.warning(
            f"Job/core number '{config['run_atacseq']['jobs']}' "
            "is set in configuration! Using..."
        )
    
    if config['run_atacseq']['profile'] is None:
        if args.profile:
            config['run_atacseq']['profile'] = os.path.abspath(args.profile)
    else:
        logger.warning(
            f"Profile '{config['run_atacseq']['profile']}' "
            "is set in configuration! Using..."
        )
    
    return config


def run_atacseq_config_check(config):
    """ Evaluates input parameters for run_atacseq module. """

    # Input directory
    if config['run_atacseq']['indir'] is None:
        logger.error("Missing fastq file input directory.")
        sys.exit(1)
    else:
        if not os.path.exists(config['run_atacseq']['indir']):
            logger.error(f"Input directory '{config['run_atacseq']['indir']}' not found.")
            sys.exit(1)
        elif not glob.glob(os.path.join(config['run_atacseq']['indir'], '*.f*q.gz')):
            logger.error(f"No fastq files found in '{config['run_atacseq']['indir']}'.")
            sys.exit(1)

    # fastp options
    blocked_opts = [
        '-i', '--in1', '-o', '--out1', '-I', '--in2', '-O', '--out2', '-m', '--merge',
        '--merged_out', '--include_unmerged', '--stdin', '--stdout', '--interleaved_in', 
        '-j', '--json', '-h', '--html', '-s', '--split', '-S', '--split_by_lines', 
        '-d', '--split_prefix_digits', '-?', '--help', '-v', '--version'
    ]
    if config['run_atacseq']['fastp_opt']:
        fastp_opts = re.split('\s+', config['run_atacseq']['fastp_opt'])
        for opt in fastp_opts:
            if opt in blocked_opts:
                logger.error(f"fastp option '{opt}' is not allowed.")
                sys.exit(1)

    # HISAT2 index and options
    if config['run_atacseq']['hisat2_ref'] is None:
        logger.error("Missing HISAT2 index.")
        sys.exit(1)
    else:
        if not os.path.exists(config['run_atacseq']['hisat2_ref'] + '.1.ht2'):
            logger.error(f"Missing HISAT2 index files at '{config['run_atacseq']['hisat2_ref']}'.")
            sys.exit(1)
    
    blocked_opts = [
        '-x', '-1', '-2', '-U', '-S', '--qseq', '-f', '-r', '-c', '--nofw', '--norc', 
        '--rna-strandness', '--summary-file', '--version', '-h', '--help'
    ]
    if config['run_atacseq']['hisat2_opt']:
        hisat2_opts = re.split('\s+', config['run_atacseq']['hisat2_opt'])
        for opt in hisat2_opts:
            if opt in blocked_opts:
                logger.error(f"HISAT2 option '{opt}' is not allowed.")
                sys.exit(1)
    
    # Exclusion list
    if config['run_atacseq']['exclusion_list']:
        if not os.path.exists(config['run_atacseq']['exclusion_list']):
            logger.error(f"Exclusion list '{config['run_atacseq']['exclusion_list']}' not found.")
            sys.exit(1)
    
    # Chromosome size file
    if config['run_atacseq']['chr_sizes'] is None:
        logger.error("Missing chromosome size file.")
        sys.exit(1)
    else:
        if not os.path.exists(config['run_atacseq']['chr_sizes']):
            logger.error(f"Chromosome size file '{config['run_atacseq']['chr_sizes']}' not found.")
            sys.exit(1)

    # Effective genome size
    if config['run_atacseq']['genome_size'] is None:
        logger.error("Missing effective genome size.")
        sys.exit(1)

    # Annotation
    if config['run_atacseq']['gtf'] is None:
        logger.error("Missing gene annotation file.")
        sys.exit(1)
    else:
        if not os.path.exists(config['run_atacseq']['gtf']):
            logger.error(f"Annotation file '{config['run_atacseq']['gtf']}' not found.")
            sys.exit(1)
    
    # Profile
    if config['run_atacseq']['profile']:
        if not os.path.exists(config['run_atacseq']['profile']):
            logger.error(f"Snakemake profile '{config['run_atacseq']['profile']}' not found.")
            sys.exit(1)

