"""Option handling and configuration validation for osync modules.

This module provides functions to set and validate configuration parameters
for each osync module (init, getref, getfq, runqc, run_rnaseq, run_atacseq).
"""

import os
import re
import glob
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from .logger import logger


# Constants - Blocked command-line options for various tools
BLOCKED_FASTP_OPTS = [
    '-i', '--in1', '-o', '--out1', '-I', '--in2', '-O', '--out2', '-m', '--merge',
    '--merged_out', '--include_unmerged', '--stdin', '--stdout', '--interleaved_in', 
    '-j', '--json', '-h', '--html', '-s', '--split', '-S', '--split_by_lines', 
    '-d', '--split_prefix_digits', '-?', '--help', '-v', '--version'
]

BLOCKED_HISAT2_OPTS = [
    '-x', '-1', '-2', '-U', '-S', '--qseq', '-f', '-r', '-c', '--nofw', '--norc', 
    '--summary-file', '--version', '-h', '--help'
]

BLOCKED_HISAT2_ATACSEQ_OPTS = BLOCKED_HISAT2_OPTS + ['--rna-strandness']

BLOCKED_KALLISTO_OPTS = [
    '-i', '--index', '-o', '--output-dir', '--plaintext', '--single'
]

# Available Kraken2 database types (October 2025)
AVAILABLE_K2DB_TYPES = [
    # Pre-built RefSeq databases
    'viral', 'minusb', 
    'standard', 'standard_8', 'standard_16',
    'pluspf', 'pluspf_8', 'pluspf_16',
    'pluspfp', 'pluspfp_8', 'pluspfp_16',
    'pracken', 'eupathdb48', 'core_nt', 'gtdb_v226',
    # Pre-built 16S databases
    '16s_greengenes_13.5', '16s_rdp_11.5', '16s_silva_132', '16s_silva_138',
    # Custom builds
    'standard_custom', 'greengenes', 'rdp', 'silva', 'nt_custom'
]


# Helper functions for common operations
def validate_file_exists(filepath: Optional[str], file_description: str) -> None:
    """Validate that a file exists, exit if not."""
    if filepath and not Path(filepath).exists():
        logger.error(f"{file_description} not found: '{filepath}'")
        sys.exit(1)


def validate_directory_exists(dirpath: Optional[str], dir_description: str) -> None:
    """Validate that a directory exists, exit if not."""
    if dirpath and not Path(dirpath).exists():
        logger.error(f"{dir_description} not found: '{dirpath}'")
        sys.exit(1)


def validate_fastq_files_exist(dirpath: str) -> None:
    """Validate that directory contains fastq files."""
    if not glob.glob(os.path.join(dirpath, '*.f*q.gz')):
        logger.error(f"No fastq files found in: '{dirpath}'")
        sys.exit(1)


def validate_command_options(options: Optional[str], blocked_opts: List[str], 
                            tool_name: str) -> None:
    """Validate that user-provided command options don't include blocked options."""
    if options:
        opts_list = re.split(r'\s+', options)
        for opt in opts_list:
            if opt in blocked_opts:
                logger.error(f"{tool_name} option '{opt}' is not allowed.")
                sys.exit(1)


def set_config_param(config: Dict[str, Any], module: str, param: str, 
                     value: Any, param_description: str) -> None:
    """Set a configuration parameter if not already set."""
    if config[module][param] is None:
        config[module][param] = value
    else:
        logger.warning(
            f"{param_description} '{config[module][param]}' "
            "is set in configuration! Using..."
        )


def abspath_if_exists(path: Optional[str]) -> Optional[str]:
    """Return absolute path if input is not None."""
    return os.path.abspath(path) if path else None


def validate_hisat2_index(index_path: str) -> None:
    """Validate that HISAT2 index files exist."""
    if not Path(f"{index_path}.1.ht2").exists():
        logger.error(f"Missing HISAT2 index files at: '{index_path}'")
        sys.exit(1)


def validate_kraken2_db(db_path: str) -> None:
    """Validate that Kraken2/Bracken database exists and has required files."""
    db_dir = Path(db_path)

    if not db_dir.exists():
        logger.error(f"Kraken2/Bracken database not found: '{db_path}'")
        sys.exit(1)

    # Check required Kraken2 files
    required_files = ['hash.k2d', 'opts.k2d', 'taxo.k2d']
    for req_file in required_files:
        if not (db_dir / req_file).exists():
            logger.error(f"Missing '{req_file}' in Kraken2 database")
            sys.exit(1)

    # Check Bracken kmer distribution files
    kmer_distribs = list(db_dir.glob('database*mers.kmer_distrib'))
    if not kmer_distribs:
        logger.error(f"No Bracken '.kmer_distrib' files found in: '{db_path}'")
        sys.exit(1)


def get_latest_ensembl_release() -> int:
    """Get the latest Ensembl release number."""
    cmd = "gget ref homo_sapiens -w dna"
    result = subprocess.check_output(
        cmd, shell=True, text=True, stderr=subprocess.DEVNULL
    )
    match = re.search(r'"ensembl_release":\s+(\d+),', result)
    return int(match.group(1))


#
# init module
#
def init_config_set(config: Dict[str, Any], args: Any) -> Dict[str, Any]:
    """Set parameters for init module."""
    module = 'init'

    set_config_param(config, module, 'hisat2_ref', 
                    abspath_if_exists(args.hisat2_ref), "HISAT2 reference file")
    set_config_param(config, module, 'hisat2_index_base', 
                    args.hisat2_index_base, "HISAT2 index basename")
    set_config_param(config, module, 'hisat2_gtf', 
                    abspath_if_exists(args.hisat2_gtf), "HISAT2 annotation file")
    set_config_param(config, module, 'rrna_ref', 
                    args.rrna_ref, "Download rRNA reference")
    set_config_param(config, module, 'exclusion_lists', 
                    args.exclusion_lists, "Download exclusion lists")
    set_config_param(config, module, 'chr_sizes_ref', 
                    abspath_if_exists(args.chr_sizes_ref), "Chromosome sizes reference file")
    set_config_param(config, module, 'kallisto_ref', 
                    abspath_if_exists(args.kallisto_ref), "Kallisto reference file")
    set_config_param(config, module, 'kallisto_index', 
                    args.kallisto_index, "Kallisto index filename")
    set_config_param(config, module, 'kraken2_db_type', 
                    args.kraken2_db_type, "Kraken2 database")

    if config[module]['bracken_read_len'] is None:
        if config[module]['kraken2_db_type']:
            config[module]['bracken_read_len'] = args.bracken_read_len
    else:
        logger.warning(
            f"Bracken read length '{config[module]['bracken_read_len']}' "
            "is set in configuration! Using..."
        )

    set_config_param(config, module, 'outdir', 
                    os.path.abspath(args.outdir), "Outdir")
    set_config_param(config, module, 'jobs', 
                    args.jobs, "Job/core number")
    set_config_param(config, module, 'profile', 
                    abspath_if_exists(args.profile), "Profile")
    set_config_param(config, module, 'latency_wait', 
                    args.latency_wait, "Latency wait")

    return config


def init_config_check(config: Dict[str, Any]) -> None:
    """Validate input parameters for init module."""
    module = 'init'

    # HISAT2 index build validation
    if config[module]['hisat2_ref']:
        validate_file_exists(config[module]['hisat2_ref'], "HISAT2 reference")

    if config[module]['hisat2_ref'] and config[module]['hisat2_index_base'] is None:
        logger.error("Missing basename for HISAT2 index")
        sys.exit(1)

    if config[module]['hisat2_gtf']:
        validate_file_exists(config[module]['hisat2_gtf'], "Annotation file")

    # Reference genome for chromosome sizes
    if config[module]['chr_sizes_ref']:
        validate_file_exists(config[module]['chr_sizes_ref'], "Chromosome sizes reference")

    # kallisto index build validation
    if config[module]['kallisto_ref']:
        validate_file_exists(config[module]['kallisto_ref'], "Kallisto reference")

    if config[module]['kallisto_ref'] and config[module]['kallisto_index'] is None:
        logger.error("Missing filename for kallisto index")
        sys.exit(1)

    # kraken2/bracken DB build validation
    if config[module]['kraken2_db_type']:
        if config[module]['kraken2_db_type'] not in AVAILABLE_K2DB_TYPES:
            logger.error(
                f"Invalid Kraken2 database type: '{config[module]['kraken2_db_type']}'.\n"
                f"Available options:\n"
                f"  RefSeq: viral, minusb, standard, standard_8, standard_16, pluspf, pluspf_8, "
                f"pluspf_16, pluspfp, pluspfp_8, pluspfp_16, pracken, eupathdb48, core_nt, gtdb_v226\n"
                f"  16S: 16s_greengenes_13.5, 16s_rdp_11.5, 16s_silva_132, 16s_silva_138\n"
                f"  Custom: standard_custom, greengenes, rdp, silva, nt_custom"
            )
            sys.exit(1)

    # Profile validation
    if config[module]['profile']:
        validate_file_exists(config[module]['profile'], "Snakemake profile")


#
# getref module
#
def getref_config_set(config: Dict[str, Any], args: Any) -> Dict[str, Any]:
    """Set parameters for getref module."""
    module = 'getref'

    if config[module]['release'] is None:
        ensembl_latest = get_latest_ensembl_release()
        config[module]['release'] = args.release if args.release is not None else ensembl_latest
    else:
        logger.warning(f"Release '{config[module]['release']}' is set in configuration! Using...")

    # List available species if requested
    if args.list_sp:
        subprocess.check_call(
            f"gget ref --list_species --release {config[module]['release']}", 
            shell=True
        )
        sys.exit(0)

    set_config_param(config, module, 'species', args.species, "Species")
    set_config_param(config, module, 'which', args.which, "Reference files")
    set_config_param(config, module, 'outdir', os.path.abspath(args.outdir), "Outdir")
    set_config_param(config, module, 'jobs', args.jobs, "Job/core number")
    set_config_param(config, module, 'profile', abspath_if_exists(args.profile), "Profile")
    set_config_param(config, module, 'latency_wait', args.latency_wait, "Latency wait")

    return config


def getref_config_check(config: Dict[str, Any]) -> None:
    """Validate input parameters for getref module."""
    module = 'getref'

    # Ensembl release validation
    ensembl_latest = get_latest_ensembl_release()
    if config[module]['release'] > ensembl_latest:
        logger.error(f"Ensembl release '{config[module]['release']}' is not valid")
        sys.exit(1)

    # Species validation
    if config[module]['species'] is None:
        logger.error("Species is not specified")
        sys.exit(1)

    # Profile validation
    if config[module]['profile']:
        validate_file_exists(config[module]['profile'], "Snakemake profile")


#
# getfq module
#
def getfq_config_set(config: Dict[str, Any], args: Any) -> Dict[str, Any]:
    """Set parameters for getfq module."""
    if 'getfq' not in config:
        config['getfq'] = {}

    # Validate mutually exclusive ID sources
    has_cmdline = hasattr(args, 'ids') and args.ids
    has_file = hasattr(args, 'ids_file') and args.ids_file

    if has_cmdline and has_file:
        logger.error("Cannot use both command-line IDs and --ids-file simultaneously")
        sys.exit(1)

    # Load IDs
    all_ids = []
    if has_file:
        if not Path(args.ids_file).exists():
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

    # Deduplicate IDs
    unique_ids = list(dict.fromkeys(all_ids))
    if len(unique_ids) < len(all_ids):
        logger.warning(f"Removed {len(all_ids) - len(unique_ids)} duplicate IDs")

    config['getfq']['ids'] = ' '.join(unique_ids)
    logger.warning(f"Total IDs: {len(unique_ids)}")

    # Set other parameters
    if hasattr(args, 'ftp'): 
        config['getfq']['ftp'] = args.ftp
    if hasattr(args, 'batch_size'): 
        config['getfq']['batch_size'] = args.batch_size
        if args.batch_size > 1: 
            logger.warning(f"Batch size: {args.batch_size}")
    if hasattr(args, 'outdir') and args.outdir: 
        config['getfq']['outdir'] = args.outdir
    if hasattr(args, 'jobs') and args.jobs: 
        config['getfq']['jobs'] = args.jobs
    if hasattr(args, 'profile'): 
        config['getfq']['profile'] = args.profile
    if hasattr(args, 'latency_wait') and args.latency_wait:
        config['getfq']['latency_wait'] = args.latency_wait
    if hasattr(args, 'dryrun'): 
        config['getfq']['dryrun'] = args.dryrun

    return config


def getfq_config_check(config: Dict[str, Any]) -> None:
    """Validate input parameters for getfq module."""
    if 'getfq' not in config:
        logger.error("getfq configuration not found")
        sys.exit(1)

    if 'ids' not in config['getfq'] or not config['getfq']['ids']:
        logger.error("No accession IDs provided")
        sys.exit(1)

    # Create output directory if needed
    if 'outdir' in config['getfq']:
        outdir = Path(config['getfq']['outdir'])
        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)
            logger.warning(f"Created output directory: {outdir}")

    # Validate batch size
    batch_size = config['getfq'].get('batch_size', 1)
    if batch_size < 1:
        logger.error(f"Invalid batch size: {batch_size}")
        sys.exit(1)

    # Profile validation
    if config['getfq'].get('profile'):
        validate_file_exists(config['getfq']['profile'], "Snakemake profile")


#
# runqc module
#
def runqc_config_set(config: Dict[str, Any], args: Any) -> Dict[str, Any]:
    """Set parameters for runqc module."""
    module = 'runqc'

    set_config_param(config, module, 'indir', abspath_if_exists(args.indir), "Input directory")

    # Handle sampling number (convert to int if > 1)
    if config[module]['number'] is None:
        number = int(args.number) if args.number > 1 else args.number
        config[module]['number'] = number
    else:
        logger.warning(f"Sampling number '{config[module]['number']}' is set in configuration! Using...")

    set_config_param(config, module, 'seed', args.seed, "Random seed")
    set_config_param(config, module, 'fastp_opt', args.fastp_opt, "fastp options")
    set_config_param(config, module, 'hisat2_ref', abspath_if_exists(args.hisat2_ref), "HISAT2 index")
    set_config_param(config, module, 'gtf', abspath_if_exists(args.gtf), "Annotation file")
    set_config_param(config, module, 'kraken2_db', abspath_if_exists(args.kraken2_db), "Kraken DB")
    set_config_param(config, module, 'outdir', os.path.abspath(args.outdir), "Outdir")
    set_config_param(config, module, 'jobs', args.jobs, "Job/core number")
    set_config_param(config, module, 'profile', abspath_if_exists(args.profile), "Profile")
    set_config_param(config, module, 'latency_wait', args.latency_wait, "Latency wait")

    return config


def runqc_config_check(config: Dict[str, Any]) -> None:
    """Validate input parameters for runqc module."""
    module = 'runqc'

    # Input directory validation
    if config[module]['indir'] is None:
        logger.error("Missing fastq file input directory")
        sys.exit(1)

    validate_directory_exists(config[module]['indir'], "Input directory")
    validate_fastq_files_exist(config[module]['indir'])

    # fastp options validation
    validate_command_options(config[module]['fastp_opt'], BLOCKED_FASTP_OPTS, "fastp")

    # HISAT2 index validation
    if config[module]['hisat2_ref']:
        validate_hisat2_index(config[module]['hisat2_ref'])

    # Annotation validation
    if config[module]['gtf']:
        validate_file_exists(config[module]['gtf'], "Annotation file")

    # kraken2/bracken DB validation
    if config[module]['kraken2_db']:
        validate_kraken2_db(config[module]['kraken2_db'])

    # Profile validation
    if config[module]['profile']:
        validate_file_exists(config[module]['profile'], "Snakemake profile")


#
# run_rnaseq module
#
def run_rnaseq_config_set(config: Dict[str, Any], args: Any) -> Dict[str, Any]:
    """Set parameters for run_rnaseq module."""
    module = 'run_rnaseq'

    set_config_param(config, module, 'indir', abspath_if_exists(args.indir), "Input directory")
    set_config_param(config, module, 'fastp_opt', args.fastp_opt, "fastp options")
    set_config_param(config, module, 'removal_ref', abspath_if_exists(args.removal_ref), 
                    "Read removal reference")
    set_config_param(config, module, 'hisat2_ref', abspath_if_exists(args.hisat2_ref), 
                    "HISAT2 index")
    set_config_param(config, module, 'hisat2_opt', args.hisat2_opt, "HISAT2 options")
    set_config_param(config, module, 'kallisto_ref', abspath_if_exists(args.kallisto_ref), 
                    "kallisto index")
    set_config_param(config, module, 'kallisto_opt', args.kallisto_opt, "kallisto options")
    set_config_param(config, module, 'gtf', abspath_if_exists(args.gtf), "Annotation file")

    if config[module]['exon_quant'] is None:
        if args.exon_quant:
            config[module]['exon_quant'] = args.exon_quant
    else:
        logger.warning(
            f"Exon-level quantification '{config[module]['exon_quant']}' "
            "is set in configuration! Using..."
        )

    set_config_param(config, module, 'genome', abspath_if_exists(args.genome), 
                    "Reference FASTA file")
    set_config_param(config, module, 'outdir', os.path.abspath(args.outdir), "Outdir")
    set_config_param(config, module, 'jobs', args.jobs, "Job/core number")
    set_config_param(config, module, 'profile', abspath_if_exists(args.profile), "Profile")
    set_config_param(config, module, 'latency_wait', args.latency_wait, "Latency wait")

    return config


def run_rnaseq_config_check(config: Dict[str, Any]) -> None:
    """Validate input parameters for run_rnaseq module."""
    module = 'run_rnaseq'

    # Input directory validation
    if config[module]['indir'] is None:
        logger.error("Missing fastq file input directory")
        sys.exit(1)

    validate_directory_exists(config[module]['indir'], "Input directory")
    validate_fastq_files_exist(config[module]['indir'])

    # fastp options validation
    validate_command_options(config[module]['fastp_opt'], BLOCKED_FASTP_OPTS, "fastp")

    # Read removal reference validation
    if config[module]['removal_ref']:
        validate_file_exists(config[module]['removal_ref'], "Read removal reference file")

    # HISAT2 index and options validation
    if config[module]['hisat2_ref'] is None:
        logger.error("Missing HISAT2 index")
        sys.exit(1)
    validate_hisat2_index(config[module]['hisat2_ref'])
    validate_command_options(config[module]['hisat2_opt'], BLOCKED_HISAT2_OPTS, "HISAT2")

    # kallisto index and options validation
    if config[module]['kallisto_ref']:
        validate_file_exists(config[module]['kallisto_ref'], "kallisto index file")
    validate_command_options(config[module]['kallisto_opt'], BLOCKED_KALLISTO_OPTS, "kallisto")

    # Annotation validation
    if config[module]['gtf'] is None:
        logger.error("Missing gene annotation file")
        sys.exit(1)
    validate_file_exists(config[module]['gtf'], "Annotation file")

    # Exon-level quantification validation
    if config[module]['exon_quant']:
        if config[module]['genome']:
            validate_file_exists(config[module]['genome'], "Reference FASTA file")
    else:
        if config[module]['genome']:
            logger.error(
                f"Reference FASTA file '{config[module]['genome']}' should be "
                "provided with -e/--exon_quant flag for exon-level quantification"
            )
            sys.exit(1)

    # Profile validation
    if config[module]['profile']:
        validate_file_exists(config[module]['profile'], "Snakemake profile")


#
# run_atacseq module
#
def run_atacseq_config_set(config: Dict[str, Any], args: Any) -> Dict[str, Any]:
    """Set parameters for run_atacseq module."""
    module = 'run_atacseq'

    set_config_param(config, module, 'indir', abspath_if_exists(args.indir), "Input directory")
    set_config_param(config, module, 'fastp_opt', args.fastp_opt, "fastp options")
    set_config_param(config, module, 'hisat2_ref', abspath_if_exists(args.hisat2_ref), 
                    "HISAT2 index")
    set_config_param(config, module, 'hisat2_opt', args.hisat2_opt, "HISAT2 options")
    set_config_param(config, module, 'exclusion_list', abspath_if_exists(args.exclusion_list), 
                    "Exclusion list")
    set_config_param(config, module, 'chr_sizes', abspath_if_exists(args.chr_sizes), 
                    "Chromosome size file")
    set_config_param(config, module, 'tn5_shift', args.tn5_shift, "Tn5 shift")
    set_config_param(config, module, 'genome_size', args.genome_size, "Effective genome size")
    set_config_param(config, module, 'gtf', abspath_if_exists(args.gtf), "Annotation file")
    set_config_param(config, module, 'outdir', os.path.abspath(args.outdir), "Outdir")
    set_config_param(config, module, 'jobs', args.jobs, "Job/core number")
    set_config_param(config, module, 'profile', abspath_if_exists(args.profile), "Profile")
    set_config_param(config, module, 'latency_wait', args.latency_wait, "Latency wait")

    return config


def run_atacseq_config_check(config: Dict[str, Any]) -> None:
    """Validate input parameters for run_atacseq module."""
    module = 'run_atacseq'

    # Input directory validation
    if config[module]['indir'] is None:
        logger.error("Missing fastq file input directory")
        sys.exit(1)

    validate_directory_exists(config[module]['indir'], "Input directory")
    validate_fastq_files_exist(config[module]['indir'])

    # fastp options validation
    validate_command_options(config[module]['fastp_opt'], BLOCKED_FASTP_OPTS, "fastp")

    # HISAT2 index and options validation
    if config[module]['hisat2_ref'] is None:
        logger.error("Missing HISAT2 index")
        sys.exit(1)
    validate_hisat2_index(config[module]['hisat2_ref'])
    validate_command_options(config[module]['hisat2_opt'], BLOCKED_HISAT2_ATACSEQ_OPTS, "HISAT2")

    # Exclusion list validation
    if config[module]['exclusion_list']:
        validate_file_exists(config[module]['exclusion_list'], "Exclusion list")

    # Chromosome size file validation
    if config[module]['chr_sizes'] is None:
        logger.error("Missing chromosome size file")
        sys.exit(1)
    validate_file_exists(config[module]['chr_sizes'], "Chromosome size file")

    # Effective genome size validation
    if config[module]['genome_size'] is None:
        logger.error("Missing effective genome size")
        sys.exit(1)

    # Annotation validation
    if config[module]['gtf'] is None:
        logger.error("Missing gene annotation file")
        sys.exit(1)
    validate_file_exists(config[module]['gtf'], "Annotation file")

    # Profile validation
    if config[module]['profile']:
        validate_file_exists(config[module]['profile'], "Snakemake profile")
