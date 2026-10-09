"""Command-line interface for omicsync.

This module provides the main entry point for the omicsync toolkit, including
argument parsing, configuration management, and workflow execution for all
supported modules.
"""

import argparse
import sys
from pathlib import Path

from .__init__ import __version__
from .logger import logger
from .config_handler import (
    make_configfile,
    load_configfile,
    update_configfile,
    clear_configfile
)
from .opt_handler import (
    init_config_set, init_config_check,
    getref_config_set, getref_config_check,
    getfq_config_set, getfq_config_check,
    runqc_config_set, runqc_config_check,
    run_rnaseq_config_set, run_rnaseq_config_check,
    run_atacseq_config_set, run_atacseq_config_check
)
from .workflow.utils.utils import (
    execute_workflow,
    handle_getref_workflow,
    handle_getfq_workflow
)


# Module registry - centralizes module configuration
MODULES = {
    'init': {
        'config_set': init_config_set,
        'config_check': init_config_check,
        'snakefile': 'init.smk',
        'description': 'Initializes resources for filtering, alignment processing and kraken2 analysis.'
    },
    'getref': {
        'config_set': getref_config_set,
        'config_check': getref_config_check,
        'snakefile': 'getref.smk',
        'description': 'Retrieves reference sequences and annotations for species.'
    },
    'getfq': {
        'config_set': getfq_config_set,
        'config_check': getfq_config_check,
        'snakefile': 'getfq.smk',
        'description': 'Retrieves fastq files from a sequence repository.'
    },
    'runqc': {
        'config_set': runqc_config_set,
        'config_check': runqc_config_check,
        'snakefile': 'runqc.smk',
        'description': 'Executes quality control module.'
    },
    'run_rnaseq': {
        'config_set': run_rnaseq_config_set,
        'config_check': run_rnaseq_config_check,
        'snakefile': 'run_rnaseq.smk',
        'description': 'Executes RNA-seq module.'
    },
    'run_atacseq': {
        'config_set': run_atacseq_config_set,
        'config_check': run_atacseq_config_check,
        'snakefile': 'run_atacseq.smk',
        'description': 'Executes ATAC-seq module.'
    }
}


def add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add common arguments to a parser.

    Args:
        parser: ArgumentParser to add arguments to
    """
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        type=str,
        help="Output directory for files (default is current directory)."
    )
    parser.add_argument(
        "-j", "--jobs",
        default=1,
        type=int,
        help="Maximum number of jobs/cores to use in parallel (default=1)."
    )
    parser.add_argument(
        "-c", "--configfile",
        type=str,
        help="Sets path to configuration file."
    )
    parser.add_argument(
        "--profile",
        type=str,
        help="Sets path to snakemake profile (useful for cluster execution)."
    )
    parser.add_argument(
        "--latency_wait",
        default=1,
        type=int,
        help="Snakemake latency wait time in seconds (default=1)."
    )
    parser.add_argument(
        "--dryrun",
        action='store_true',
        default=None,
        help="Performs a dry run without execution."
    )


def create_argument_parser() -> argparse.ArgumentParser:
    """Create and configure the main argument parser.

    Returns:
        Configured ArgumentParser instance
    """
    parser = argparse.ArgumentParser(
        description=(
            f"omicsync v{__version__}: A versatile toolkit for efficient processing "
            "and integration of omics data."
        ),
        add_help=False
    )

    parser.add_argument("-v", "--version", action="store_true", help="Prints version.")
    parser.add_argument("-h", "--help", action="store_true", help="Prints help message.")

    return parser


def setup_subparsers(parser: argparse.ArgumentParser) -> dict:
    """Set up all subparsers for omicsync commands.

    Args:
        parser: Main ArgumentParser instance

    Returns:
        Dictionary mapping command names to their parsers
    """
    parent = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest="command")
    parsers = {}

    # makeconfig subparser
    msg = "Makes configuration file for omicsync execution."
    parsers['makeconfig'] = subparsers.add_parser(
        "makeconfig", parents=[parent], description=msg, help=msg, add_help=True
    )
    parsers['makeconfig'].add_argument(
        "-o", "--output", type=str, default="config.yaml",
        help="Indicates output file for configuration (default='./config.yaml')."
    )

    # clearconfig subparser
    msg = "Empties parameters set in the configuration file."
    parsers['clearconfig'] = subparsers.add_parser(
        "clearconfig", parents=[parent], description=msg, help=msg, add_help=True
    )
    parsers['clearconfig'].add_argument(
        "-w", "--which", type=str,
        help=(
            "Indicates which module to clear in the configuration file. "
            "Available options: init, getref, getfq, runqc, run_rnaseq, run_atacseq, all."
        )
    )
    parsers['clearconfig'].add_argument(
        "-c", "--configfile", type=str, help="Sets path to configuration file."
    )

    # init subparser
    parsers['init'] = subparsers.add_parser(
        "init", parents=[parent], 
        description=MODULES['init']['description'],
        help=MODULES['init']['description'],
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parsers['init'].add_argument(
        "--hisat2_ref", type=str, 
        help="Path to reference sequence file for HISAT2 index."
    )
    parsers['init'].add_argument(
        "--hisat2_index_base", type=str, 
        help="Basename to use for HISAT2 index files."
    )
    parsers['init'].add_argument(
        "--hisat2_gtf", type=str,
        help=(
            "Path to reference annotation file. Includes exon and splice site information "
            "for RNA-seq analysis. Requires ~160Gb RAM for human genome."
        )
    )
    parsers['init'].add_argument(
        "--rrna_ref", action='store_true', default=None, 
        help="Downloads rRNA reference files to output directory."
    )
    parsers['init'].add_argument(
        "--exclusion_lists", action='store_true', default=None, 
        help="Downloads exclusion list files to output directory."
    )
    parsers['init'].add_argument(
        "--chr_sizes_ref", type=str, 
        help="Path to reference sequence file for retrieving chromosome sizes."
    )
    parsers['init'].add_argument(
        "--kallisto_ref", type=str, 
        help="Path to reference sequence file for kallisto index."
    )
    parsers['init'].add_argument(
        "--kallisto_index", type=str, 
        help="Filename for kallisto index."
    )
    parsers['init'].add_argument(
        "--kraken2_db_type", type=str,
        metavar='DB_TYPE',
        help=(
            "Type of Kraken2 database. "
            "PRE-BUILT REFSEQ (Oct 2025, includes Bracken 50/75/100/150/200/250/300-mers): "
            "viral (0.6GB), minusb (11GB), standard (94GB), standard_8 (7.5GB), "
            "standard_16 (14.9GB), pluspf (100.6GB), pluspf_8 (7.5GB), pluspf_16 (14.9GB), "
            "pluspfp (214.5GB), pluspfp_8 (7.5GB), pluspfp_16 (14.9GB), pracken (351.6GB), "
            "eupathdb48 (11GB), core_nt (316.2GB), gtdb_v226 (644GB). "
            "PRE-BUILT 16S (Mar 2020, includes Bracken 100/150/200-mers): "
            "16s_greengenes_13.5 (73MB), 16s_rdp_11.5 (168MB), "
            "16s_silva_132 (117MB), 16s_silva_138 (112MB). "
            "CUSTOM BUILDS (requires --bracken_read_len): standard_custom, greengenes, rdp, silva, nt_custom. "
            "Capped DBs (-8/-16) trade accuracy for reduced memory. "
            "Credit: https://benlangmead.github.io/aws-indexes/k2"
        )
    )
    parsers['init'].add_argument(
        "--bracken_read_len", default=100, type=int, 
        help="Read length for Bracken custom builds (default=100)."
    )
    add_common_args(parsers['init'])

    # getref subparser
    parsers['getref'] = subparsers.add_parser(
        "getref", parents=[parent],
        description=MODULES['getref']['description'],
        help=MODULES['getref']['description'],
        add_help=True
    )
    parsers['getref'].add_argument(
        "-s", "--species", type=str, 
        help="Indicates species for file retrieval (display available with `--list_sp`)."
    )
    parsers['getref'].add_argument(
        "-l", "--list_sp", action="store_true", default=False, 
        help="Lists available species from Ensembl (specify release with `--release`)."
    )
    parsers['getref'].add_argument(
        "-w", "--which", type=str, default="dna,gtf", 
        help="Files to fetch with `gget` (default='dna,gtf'). Options: gtf, dna, cdna, ncrna, cds, pep, all."
    )
    parsers['getref'].add_argument(
        "-r", "--release", default=None, type=int, 
        help="Ensembl release to retrieve from (default=latest)."
    )
    add_common_args(parsers['getref'])

    # getfq subparser
    parsers['getfq'] = subparsers.add_parser(
        "getfq", parents=[parent],
        description=MODULES['getfq']['description'],
        help=MODULES['getfq']['description'],
        add_help=True
    )
    parsers['getfq'].add_argument(
        "ids", type=str, nargs="*", 
        help="One or multiple accession IDs for file retrieval."
    )
    parsers['getfq'].add_argument(
        "-i", "--ids_file", type=str, 
        help="Path to text file containing accession IDs (one per line). Use for large projects."
    )
    parsers['getfq'].add_argument(
        "--ftp", default=False, action='store_true', 
        help="Limits metadata retrieval to FTP download links."
    )
    parsers['getfq'].add_argument(
        "--batch_size", type=int, default=1, 
        help="Number of accession IDs to process per ffq call. Default=1 (all at once). Use >1 to split into batches."
    )
    add_common_args(parsers['getfq'])

    # runqc subparser
    parsers['runqc'] = subparsers.add_parser(
        "runqc", parents=[parent],
        description=MODULES['runqc']['description'],
        help=MODULES['runqc']['description'],
        add_help=True
    )
    parsers['runqc'].add_argument(
        "-i", "--indir", type=str, 
        help="Input directory containing fastq files."
    )
    parsers['runqc'].add_argument(
        "-n", "--number", default=200000, type=float, 
        help="Random sampling of fastq files. Can be a fraction (0-1), all reads (1), or specific number (>1) (default=200000)."
    )
    parsers['runqc'].add_argument(
        "-s", "--seed", default=23, type=int, 
        help="Random seed for sampling (default=23)."
    )
    parsers['runqc'].add_argument(
        "-f", "--fastp_opt", default="", type=str, 
        help="Additional fastp options. Do not include -i, -o, -I, -O, -j, -h."
    )
    parsers['runqc'].add_argument(
        "-x", "--hisat2_ref", default=None, type=str, 
        help="Path and basename of HISAT2 index files."
    )
    parsers['runqc'].add_argument(
        "-g", "--gtf", default=None, type=str, 
        help="Path to annotation file in gtf format."
    )
    parsers['runqc'].add_argument(
        "-k", "--kraken2_db", default=None, type=str, 
        help="Directory containing kraken2 database."
    )
    add_common_args(parsers['runqc'])

    # run_rnaseq subparser
    parsers['run_rnaseq'] = subparsers.add_parser(
        "run_rnaseq", parents=[parent],
        description=MODULES['run_rnaseq']['description'],
        help=MODULES['run_rnaseq']['description'],
        add_help=True
    )
    parsers['run_rnaseq'].add_argument(
        "-i", "--indir", type=str, 
        help="Input directory containing fastq files."
    )
    parsers['run_rnaseq'].add_argument(
        "-f", "--fastp_opt", default="", type=str, 
        help="Additional fastp options. Do not include -i, -o, -I, -O, -j, -h."
    )
    parsers['run_rnaseq'].add_argument(
        "-r", "--removal_ref", default=None, type=str, 
        help="Path to reference file for read removal."
    )
    parsers['run_rnaseq'].add_argument(
        "-x", "--hisat2_ref", default=None, type=str, 
        help="Path and basename of HISAT2 index files."
    )
    parsers['run_rnaseq'].add_argument(
        "-a", "--hisat2_opt", default="", type=str, 
        help="Additional HISAT2 options. Do not include -x, -1, -2, -U, -S, --summary-file. Provide strandness with --rna-strandness (default=unstranded)."
    )
    parsers['run_rnaseq'].add_argument(
        "-t", "--kallisto_ref", default=None, type=str, 
        help="Filename of kallisto index file."
    )
    parsers['run_rnaseq'].add_argument(
        "-k", "--kallisto_opt", default="", type=str, 
        help="Additional kallisto options. Do not include -i, -o, --single. Provide strandness with --fr-stranded or --rf-stranded."
    )
    parsers['run_rnaseq'].add_argument(
        "-g", "--gtf", default=None, type=str, 
        help="Path to annotation file in gtf format."
    )
    parsers['run_rnaseq'].add_argument(
        "-e", "--exon_quant", action='store_true', default=None, 
        help="Performs exon-level quantification with featureCounts."
    )
    parsers['run_rnaseq'].add_argument(
        "-G", "--genome", default=None, type=str, 
        help="Path to FASTA-format file with reference sequences used in read mapping."
    )
    add_common_args(parsers['run_rnaseq'])

    # run_atacseq subparser
    parsers['run_atacseq'] = subparsers.add_parser(
        "run_atacseq", parents=[parent],
        description=MODULES['run_atacseq']['description'],
        help=MODULES['run_atacseq']['description'],
        add_help=True
    )
    parsers['run_atacseq'].add_argument(
        "-i", "--indir", type=str, 
        help="Input directory containing fastq files."
    )
    parsers['run_atacseq'].add_argument(
        "-f", "--fastp_opt", default="", type=str, 
        help="Additional fastp options. Do not include -i, -o, -I, -O, -j, -h."
    )
    parsers['run_atacseq'].add_argument(
        "-x", "--hisat2_ref", default=None, type=str, 
        help="Path and basename of HISAT2 index files."
    )
    parsers['run_atacseq'].add_argument(
        "-a", "--hisat2_opt", default="", type=str, 
        help="Additional HISAT2 options. Do not include -x, -1, -2, -U, -S, --summary-file."
    )
    parsers['run_atacseq'].add_argument(
        "-e", "--exclusion_list", default=None, type=str, 
        help="Path to genome exclusion list in bed format."
    )
    parsers['run_atacseq'].add_argument(
        "-l", "--chr_sizes", default=None, type=str, 
        help="Path to tab-delimited file with chromosome sizes."
    )
    parsers['run_atacseq'].add_argument(
        "-t", "--tn5_shift", action='store_true', default=None, 
        help="Adjusts read mapping for Tn5 transposase shift."
    )
    parsers['run_atacseq'].add_argument(
        "-m", "--genome_size", default=None, type=int, 
        help="Effective genome size. Reference values at https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html"
    )
    parsers['run_atacseq'].add_argument(
        "-g", "--gtf", default=None, type=str, 
        help="Path to annotation file in gtf format."
    )
    add_common_args(parsers['run_atacseq'])

    return parsers


def process_module_workflow(module: str, args) -> None:
    """Process a standard module workflow.

    Args:
        module: Module name
        args: Parsed command-line arguments
    """
    config = load_configfile(args.configfile)
    config = MODULES[module]['config_set'](config, args)
    MODULES[module]['config_check'](config)
    update_configfile(args.configfile, config)

    # Special handling for modules with unique requirements
    if module == 'getref':
        handle_getref_workflow(config, args)
    elif module == 'getfq':
        handle_getfq_workflow(config, args)
    else:
        execute_workflow(module, config, args.configfile, args.dryrun)


def main() -> None:
    """Main entry point for the omicsync command-line interface."""
    parser = create_argument_parser()
    subparser_dict = setup_subparsers(parser)
    args = parser.parse_args()

    # Display help for main parser
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and args.help):
        parser.print_help(sys.stderr)
        sys.exit(0)

    # Version return
    if args.version:
        print(f"omicsync v{__version__}")
        sys.exit(0)

    # Display subparser specific help when only subparser command is specified
    if len(sys.argv) == 2 and sys.argv[1] in subparser_dict:
        if args.command != "makeconfig":
            subparser_dict[sys.argv[1]].print_help(sys.stderr)
            sys.exit(0)

    # Handle commands
    if args.command == "makeconfig":
        make_configfile(args.output)

    elif args.command == "clearconfig":
        clear_configfile(args.configfile, args.which)

    elif args.command in MODULES:
        process_module_workflow(args.command, args)

    else:
        parser.print_help(sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
