import argparse
import os
import sys
import json
import subprocess

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


def main():
    """
    Entry point for command-line interface and argument parsers. 
    """

    #
    # Main parser
    #
    parser = argparse.ArgumentParser(
        description=(
            f"fqseek v{__version__}: an integrated framework for reproducible "
            "and efficient processing of high-throughput sequencing data."
        ),
        add_help=False
    )

    # Fix for parsing error
    parent = argparse.ArgumentParser(add_help=False)

    # Add version
    parser.add_argument(
        "-v", 
        "--version", 
        action="store_true", 
        help="Prints version."
    )
    # Add help
    parser.add_argument(
        "-h", 
        "--help", 
        action="store_true", 
        help="Prints help message."
    )

    #
    # Subparsers
    #
    subparsers = parser.add_subparsers(dest="command")

    # makeconfig subparser
    makeconfig_msg = "Makes configuration file for fqseek execution."
    makeconfig_parser = subparsers.add_parser(
        "makeconfig", 
        parents=[parent], 
        description=makeconfig_msg, 
        help=makeconfig_msg, 
        add_help=True
    )

    makeconfig_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="config.yaml",
        required=False,
        help="Indicates output file for configuration (default='./config.yaml')."
    )


    # clearconfig subparser
    clearconfig_msg = "Empties parameters set in the configuration file."
    clearconfig_parser = subparsers.add_parser(
        "clearconfig", 
        parents=[parent], 
        description=clearconfig_msg, 
        help=clearconfig_msg, 
        add_help=True
    )

    clearconfig_parser.add_argument(
        "-w",
        "--which",
        type=str,
        required=False,
        help=(
            "Indicates which module to clear in the configuration file. "
            "Available options are: init, getref, getfq, runqc, run_rnaseq, "
            "run_atacseq, all."
        )
    )
    clearconfig_parser.add_argument(
        "-c",
        "--configfile",
        type=str,
        required=False,
        help="Sets path to configuration file."
    )


    # init subparser
    init_msg = "Initializes resources for filtering, alignment processing and kraken2 analysis."
    init_parser = subparsers.add_parser(
        "init", 
        parents=[parent], 
        description=init_msg, 
        help=init_msg, 
        add_help=True
    )

    init_parser.add_argument(
        "--hisat2_ref",
        type=str,
        required=False,
        help="Path to reference sequence file for HISAT2 index."
    )
    init_parser.add_argument(
        "--hisat2_index_base",
        type=str,
        required=False,
        help="Basename to use for HISAT2 index files."
    )
    init_parser.add_argument(
        "--hisat2_gtf",
        type=str,
        required=False,
        help=(
            "Path to reference annotation file. This includes exon and splice site "
            "information for RNA-seq analysis. Requires ~160Gb RAM for human genome."
        )
    )
    init_parser.add_argument(
        "--rrna_ref",
        action='store_true',
        default=None,
        required=False,
        help="Downloads rRNA reference files to output directory."
    )
    init_parser.add_argument(
        "--exclusion_lists",
        action='store_true',
        default=None,
        required=False,
        help="Downloads exclusion list files to output directory."
    )
    init_parser.add_argument(
        "--chr_sizes_ref",
        type=str,
        required=False,
        help="Path to reference sequence file for retrieving chromosome sizes."
    )
    init_parser.add_argument(
        "--kallisto_ref",
        type=str,
        required=False,
        help="Path to reference sequence file for kallisto index."
    )
    init_parser.add_argument(
        "--kallisto_index",
        type=str,
        required=False,
        help="Filename for kallisto index."
    )
    init_parser.add_argument(
        "--kraken2_db_type",
        type=str,
        required=False,
        help=(
            "Type of kraken2 database to build. Choose from: 'standard', 'greengenes', "
            "'rdp', 'silva', 'nt', 'std16_prebuilt' and 'plusPFP16_prebuilt'. Make sure "
            "enough RAM and disk space are available for 'standard' and 'nt' build. Additional "
            "prebuilt databases available at https://benlangmead.github.io/aws-indexes/k2."
        )
    )
    init_parser.add_argument(
        "--bracken_read_len",
        default=100,
        type=int,
        required=False,
        help="Read length for classifications (default=100)."
    )
    init_parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        type=str,
        required=False,
        help="Output directory for files (default is current directory)."
    )
    init_parser.add_argument(
        "-j",
        "--jobs",
        default=1,
        type=int,
        required=False,
        help="Maximum number of jobs/cores to use in parallel (default=1)."
    )
    init_parser.add_argument(
        "-c",
        "--configfile",
        type=str,
        required=False,
        help="Sets path to configuration file."
    )
    init_parser.add_argument(
        "-p",
        "--profile",
        type=str,
        required=False,
        help="Sets path to snakemake profile (useful for cluster execution)."
    )
    init_parser.add_argument(
        "-z",
        "--dryrun",
        action='store_true',
        default=None,
        required=False,
        help="Performs a dry run without execution."
    )


    # getref subparser
    getref_msg = "Retrieves reference sequences and annotations for species."
    getref_parser = subparsers.add_parser(
        "getref", 
        parents=[parent], 
        description=getref_msg, 
        help=getref_msg, 
        add_help=True
    )

    getref_parser.add_argument(
        "-s",
        "--species",
        type=str,
        default=None,
        required=False,
        help="Indicates species for file retrieval (display available with `--list_sp`)."
    )
    getref_parser.add_argument(
        "-l",
        "--list_sp",
        action="store_true",
        default=False,
        required=False,
        help="Lists available species from Ensembl (specify release with `--release`)."
    )
    getref_parser.add_argument(
        "-w",
        "--which",
        type=str,
        default="dna,gtf",
        required=False,
        help=(
            "Files to fetch with `gget` (default='dna,gtf'). Available "
            "options are: gtf, dna, cdna, ncrna, cds, pep, all."
        )
    )
    getref_parser.add_argument(
        "-r",
        "--release",
        default=None,
        type=int,
        required=False,
        help="Ensembl release to retrieve from (default=latest)."
    )
    getref_parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        type=str,
        required=False,
        help="Output directory for files (default is current directory)."
    )
    getref_parser.add_argument(
        "-j",
        "--jobs",
        default=1,
        type=int,
        required=False,
        help="Maximum number of jobs/cores to use in parallel (default=1)."
    )
    getref_parser.add_argument(
        "-c",
        "--configfile",
        type=str,
        required=False,
        help="Sets path to configuration file."
    )
    getref_parser.add_argument(
        "-p",
        "--profile",
        type=str,
        required=False,
        help="Sets path to snakemake profile (useful for cluster execution)."
    )
    getref_parser.add_argument(
        "-z",
        "--dryrun",
        action='store_true',
        default=None,
        required=False,
        help="Performs a dry run without execution."
    )


    # getfq subparser
    getfq_msg = "Retrieves fastq files from a sequence respository."
    getfq_parser = subparsers.add_parser(
        "getfq", 
        parents=[parent], 
        description=getfq_msg, 
        help=getfq_msg, 
        add_help=True
    )

    getfq_parser.add_argument(
        "ids",
        type=str,
        nargs="*",
        help="One or multiple accession IDs for file retrieval."
    )
    getfq_parser.add_argument(
        "-i",
        "--ids_file",
        type=str,
        required=False,
        help=(
            "Path to text file containing accession IDs (one per line). "
            "Use this for large projects with many IDs instead of command-line input."
        )
    )
    getfq_parser.add_argument(
        "--ftp",
        default=False,
        action='store_true',
        help="Limits metadata retrieval to FTP download links."
    )
    getfq_parser.add_argument(
        "--batch_size",
        type=int,
        default=1,
        required=False,
        help=(
            "Number of accession IDs to process per ffq call. "
            "Default=1 (all at once). Use >1 to split into batches."
        )
    )
    getfq_parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        type=str,
        required=False,
        help="Output directory for files (default is current directory)."
    )
    getfq_parser.add_argument(
        "-j",
        "--jobs",
        default=1,
        type=int,
        required=False,
        help="Maximum number of jobs/cores to use in parallel (default=1)."
    )
    getfq_parser.add_argument(
        "-c",
        "--configfile",
        type=str,
        required=False,
        help="Sets path to configuration file."
    )
    getfq_parser.add_argument(
        "-p",
        "--profile",
        type=str,
        required=False,
        help="Sets path to snakemake profile (useful for cluster execution)."
    )
    getfq_parser.add_argument(
        "-z",
        "--dryrun",
        action='store_true',
        default=None,
        required=False,
        help="Performs a dry run without execution."
    )


    # runqc subparser
    runqc_msg = "Executes quality control module."
    runqc_parser = subparsers.add_parser(
        "runqc", 
        parents=[parent], 
        description=runqc_msg, 
        help=runqc_msg, 
        add_help=True
    )

    runqc_parser.add_argument(
        "-i",
        "--indir",
        type=str,
        required=False,
        help="Input directory containing fastq files."
    )
    runqc_parser.add_argument(
        "-n",
        "--number",
        default=200000,
        type=float,
        required=False,
        help=(
            "Specifies random sampling of fastq files. Can be a fraction of reads (0-1), "
            "all reads (1) or a specific number of reads (>1) (default=200000)."
        )
    )
    runqc_parser.add_argument(
        "-s",
        "--seed",
        default=23,
        type=int,
        required=False,
        help="Specifies random seed used for sampling (default=23)."
    )
    runqc_parser.add_argument(
        "-f",
        "--fastp_opt",
        default="",
        type=str,
        required=False,
        help=(
            "Additional options for fastq file processing with fastp. Do not include "
            "arguments for input and output files (-i, -o, -I, -O, -j, -h)."
        )
    )
    runqc_parser.add_argument(
        "-x",
        "--hisat2_ref",
        default=None,
        type=str,
        required=False,
        help="Path and basename of HISAT2 index files."
    )
    runqc_parser.add_argument(
        "-g",
        "--gtf",
        default=None,
        type=str,
        required=False,
        help="Path to annotation file in gtf format."
    )
    runqc_parser.add_argument(
        "-k",
        "--kraken2_db",
        default=None,
        type=str,
        required=False,
        help="Directory containing database for kraken2 analysis."
    )
    runqc_parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        type=str,
        required=False,
        help="Output directory for files (default is current directory)."
    )
    runqc_parser.add_argument(
        "-j",
        "--jobs",
        default=1,
        type=int,
        required=False,
        help="Maximum number of jobs/cores to use in parallel (default=1)."
    )
    runqc_parser.add_argument(
        "-c",
        "--configfile",
        type=str,
        required=False,
        help="Sets path to configuration file."
    )
    runqc_parser.add_argument(
        "-p",
        "--profile",
        type=str,
        required=False,
        help="Sets path to snakemake profile (useful for cluster execution)."
    )
    runqc_parser.add_argument(
        "-z",
        "--dryrun",
        action='store_true',
        default=None,
        required=False,
        help="Performs a dry run without execution."
    )

    # run_rnaseq subparser
    run_rnaseq_msg = "Executes RNA-seq module."
    run_rnaseq_parser = subparsers.add_parser(
        "run_rnaseq", 
        parents=[parent], 
        description=run_rnaseq_msg, 
        help=run_rnaseq_msg, 
        add_help=True
    )

    run_rnaseq_parser.add_argument(
        "-i",
        "--indir",
        type=str,
        required=False,
        help="Input directory containing fastq files."
    )
    run_rnaseq_parser.add_argument(
        "-f",
        "--fastp_opt",
        default="",
        type=str,
        required=False,
        help=(
            "Additional options for fastq file processing with fastp. Do not include "
            "arguments: -i, -o, -I, -O, -j, -h."
        )
    )
    run_rnaseq_parser.add_argument(
        "-r",
        "--removal_ref",
        default=None,
        type=str,
        required=False,
        help="Path to reference file for read removal."
    )
    run_rnaseq_parser.add_argument(
        "-x",
        "--hisat2_ref",
        default=None,
        type=str,
        required=False,
        help="Path and basename of HISAT2 index files."
    )
    run_rnaseq_parser.add_argument(
        "-a",
        "--hisat2_opt",
        default="",
        type=str,
        required=False,
        help=(
            "Additional options for read alignment with HISAT2. Do not include "
            "arguments: -x, -1, -2, -U, -S, --summary-file. Provide strandness with "
            "--rna-strandness (default=unstranded)."
        )
    )
    run_rnaseq_parser.add_argument(
        "-t",
        "--kallisto_ref",
        default=None,
        type=str,
        required=False,
        help="Filename of kallisto index file."
    )
    run_rnaseq_parser.add_argument(
        "-k",
        "--kallisto_opt",
        default="",
        type=str,
        required=False,
        help=(
            "Additional options for quantification with kallisto. Do not include "
            "arguments: -i, -o, --single. Provide strandness with --fr-stranded "
            "or --rf-stranded. For single-end provide estimated fragment length (-l) "
            "and fragment length standard deviation (-s)."
        )
    )
    run_rnaseq_parser.add_argument(
        "-g",
        "--gtf",
        default=None,
        type=str,
        required=False,
        help="Path to annotation file in gtf format."
    )
    run_rnaseq_parser.add_argument(
        "-e",
        "--exon_quant",
        action='store_true',
        default=None,
        required=False,
        help="Performs exon-level quantification with featureCounts."
    )
    run_rnaseq_parser.add_argument(
        "-G",
        "--genome",
        default=None,
        type=str,
        required=False,
        help="Path to FASTA-format file with reference sequences used in read mapping."
    )
    run_rnaseq_parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        type=str,
        required=False,
        help="Output directory for files (default is current directory)."
    )
    run_rnaseq_parser.add_argument(
        "-j",
        "--jobs",
        default=1,
        type=int,
        required=False,
        help="Maximum number of jobs/cores to use in parallel (default=1)."
    )
    run_rnaseq_parser.add_argument(
        "-c",
        "--configfile",
        type=str,
        required=False,
        help="Sets path to configuration file."
    )
    run_rnaseq_parser.add_argument(
        "-p",
        "--profile",
        type=str,
        required=False,
        help="Sets path to snakemake profile (useful for cluster execution)."
    )
    run_rnaseq_parser.add_argument(
        "-z",
        "--dryrun",
        action='store_true',
        default=None,
        required=False,
        help="Performs a dry run without execution."
    )

    # run_atacseq subparser
    run_atacseq_msg = "Executes ATAC-seq module."
    run_atacseq_parser = subparsers.add_parser(
        "run_atacseq", 
        parents=[parent], 
        description=run_atacseq_msg, 
        help=run_atacseq_msg, 
        add_help=True
    )

    run_atacseq_parser.add_argument(
        "-i",
        "--indir",
        type=str,
        required=False,
        help="Input directory containing fastq files."
    )
    run_atacseq_parser.add_argument(
        "-f",
        "--fastp_opt",
        default="",
        type=str,
        required=False,
        help=(
            "Additional options for fastq file processing with fastp. Do not include "
            "arguments: -i, -o, -I, -O, -j, -h."
        )
    )
    run_atacseq_parser.add_argument(
        "-x",
        "--hisat2_ref",
        default=None,
        type=str,
        required=False,
        help="Path and basename of HISAT2 index files."
    )
    run_atacseq_parser.add_argument(
        "-a",
        "--hisat2_opt",
        default="",
        type=str,
        required=False,
        help=(
            "Additional options for read alignment with HISAT2. Do not include "
            "arguments: -x, -1, -2, -U, -S, --summary-file."
        )
    )
    run_atacseq_parser.add_argument(
        "-e",
        "--exclusion_list",
        default=None,
        type=str,
        required=False,
        help="Path to genome exclusion list in bed format."
    )
    run_atacseq_parser.add_argument(
        "-l",
        "--chr_sizes",
        default=None,
        type=str,
        required=False,
        help="Path to tab-delimited file with chromosome sizes."
    )
    run_atacseq_parser.add_argument(
        "-t",
        "--tn5_shift",
        action='store_true',
        default=None,
        required=False,
        help="Adjusts read mapping for Tn5 transposase shift."
    )
    run_atacseq_parser.add_argument(
        "-m",
        "--genome_size",
        default=None,
        type=int,
        required=False,
        help=(
            "Effective genome size. Reference values available at "
            "https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html."
        )
    )
    run_atacseq_parser.add_argument(
        "-g",
        "--gtf",
        default=None,
        type=str,
        required=False,
        help="Path to annotation file in gtf format."
    )
    run_atacseq_parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        type=str,
        required=False,
        help="Output directory for files (default is current directory)."
    )
    run_atacseq_parser.add_argument(
        "-j",
        "--jobs",
        default=1,
        type=int,
        required=False,
        help="Maximum number of jobs/cores to use in parallel (default=1)."
    )
    run_atacseq_parser.add_argument(
        "-c",
        "--configfile",
        type=str,
        required=False,
        help="Sets path to configuration file."
    )
    run_atacseq_parser.add_argument(
        "-p",
        "--profile",
        type=str,
        required=False,
        help="Sets path to snakemake profile (useful for cluster execution)."
    )
    run_atacseq_parser.add_argument(
        "-z",
        "--dryrun",
        action='store_true',
        default=None,
        required=False,
        help="Performs a dry run without execution."
    )


    args = parser.parse_args()

    # Display help for main parser
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and args.help):
        parser.print_help(sys.stderr)
        sys.exit(0)
    
    # Version return
    if args.version:
        print(f"fqseek v{__version__}")
        sys.exit(0)
    
    # Display subparser specific help messages when only subparser command is specified
    subparser_cmd_dict = {
        "makeconfig" : makeconfig_parser,
        "clearconfig" : clearconfig_parser,
        "init" : init_parser,
        "getref" : getref_parser,
        "getfq" : getfq_parser,
        "runqc" : runqc_parser,
        "run_rnaseq" : run_rnaseq_parser,
        "run_atacseq" : run_atacseq_parser
    }
    if len(sys.argv) == 2 and sys.argv[1] in subparser_cmd_dict:
        if args.command != "makeconfig":
            subparser_cmd_dict[sys.argv[1]].print_help(sys.stderr)
            sys.exit(0)


    if args.command == "makeconfig":
        make_configfile(args.output)
    

    if args.command == "clearconfig":
        clear_configfile(args.configfile, args.which)
    

    if args.command == "init":
        
        # Check inputs and update configuration
        config = load_configfile(args.configfile)
        config = init_config_set(config, args)
        init_config_check(config)
        update_configfile(args.configfile, config)

        # Run init workflow
        profile = config['init']['profile']
        cmd = (
            " snakemake {profile} {dryrun} {jobs} {cores} --snakefile {snakefile} "
            " --directory {directory} --configfile {configfile} --nolock --latency-wait 60 "
            " --rerun-incomplete --rerun-triggers mtime --scheduler greedy "
        ).format(
            profile = f"--profile {profile}" if profile is not None else "",
            dryrun = f"--dry-run" if args.dryrun else "",
            jobs = f"--jobs {config['init']['jobs']}",
            cores = f"--cores {config['init']['jobs']}",
            snakefile = os.path.dirname(os.path.abspath(__file__))+"/workflow/rules/init.smk",
            directory = os.path.abspath(config['init']['outdir']),
            configfile = os.path.abspath(args.configfile)
        )

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            logger.error(e)
            sys.exit(1)
        

    if args.command == "getref":
        
        # Check inputs and update configuration
        config = load_configfile(args.configfile)
        config = getref_config_set(config, args)
        getref_config_check(config)
        update_configfile(args.configfile, config)

        # Initialize directory for download
        json_path = (
            "{outdir}/{species}_ensembl_{release}.json"
        ).format(
            species = config['getref']['species'],
            release = config['getref']['release'],
            outdir = config['getref']['outdir']
        )
        if os.path.exists(json_path):
            logger.warning(f"File exists '{json_path}'! Overwriting...")

        cmd = (
            " gget ref {species} --which {which} --release {release} "
            " --out {outdir}/{species}_ensembl_{release}.json "
        ).format(
            species = config['getref']['species'],
            which = config['getref']['which'],
            release = config['getref']['release'],
            outdir = config['getref']['outdir']
        )
        subprocess.check_call(cmd, shell=True, stderr=subprocess.DEVNULL)

        # Run getref workflow
        profile = config['getref']['profile']
        cmd = (
            " snakemake {profile} {dryrun} {jobs} {cores} --snakefile {snakefile} "
            " --directory {directory} --configfile {configfile} --nolock --latency-wait 60 "
            " --rerun-incomplete --rerun-triggers mtime --scheduler greedy "
        ).format(
            profile = f"--profile {profile}" if profile is not None else "",
            dryrun = f"--dry-run" if args.dryrun else "",
            jobs = f"--jobs {config['getref']['jobs']}",
            cores = f"--cores {config['getref']['jobs']}",
            snakefile = os.path.dirname(os.path.abspath(__file__))+"/workflow/rules/getref.smk",
            directory = os.path.abspath(config['getref']['outdir']),
            configfile = os.path.abspath(args.configfile)
        )

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            logger.error(e)
            sys.exit(1)


    if args.command == "getfq":
        
        # Check inputs and update configuration
        config = load_configfile(args.configfile)
        config = getfq_config_set(config, args)
        getfq_config_check(config)
        update_configfile(args.configfile, config)

        # Initialize directory for download
        json_path = f"{config['getfq']['outdir']}/download_meta.json"
        batch_size = config['getfq'].get('batch_size', 1)
        
        if os.path.exists(json_path):
            logger.warning(f"File exists '{json_path}'! Using...")
        else:
            ids_list = config['getfq']['ids'].split()
            total_ids = len(ids_list)
            ftp_flag = '--ftp' if config['getfq'].get('ftp', False) else ''
            
            logger.warning(f"Fetching metadata for {total_ids} accession(s)")
            
            if batch_size == 1:
                # Single ffq call for all IDs
                cmd = f"ffq {config['getfq']['ids']} {ftp_flag} -o {json_path}"
                subprocess.check_call(
                    f"""for i in {{1..5}}; do
                        echo "Retrieving download links (attempt $i/5)"
                        {cmd}
                        sleep 5 && [[ -s {json_path} ]] && exit 0 || sleep 5
                    done || (echo "ffq failed all retries" && exit 1)""",
                    shell=True, stderr=subprocess.DEVNULL
                )
            else:
                # Batch processing
                all_metadata = []
                num_batches = (total_ids + batch_size - 1) // batch_size
                
                for i in range(0, total_ids, batch_size):
                    batch = ids_list[i:i+batch_size]
                    batch_num = i // batch_size + 1
                    temp_json = f"{config['getfq']['outdir']}/batch_{batch_num}_temp.json"
                    
                    logger.warning(f"Batch {batch_num}/{num_batches}: {len(batch)} ID(s)")
                    cmd = f"ffq {' '.join(batch)} {ftp_flag} -o {temp_json}"
                    
                    try:
                        subprocess.check_call(
                            f"""for i in {{1..5}}; do
                                echo "  Retrieving download links (attempt $i/5)"
                                {cmd}
                                sleep 5 && [[ -s {temp_json} ]] && exit 0 || sleep 5
                            done || exit 1""",
                            shell=True, stderr=subprocess.DEVNULL
                        )
                        
                        with open(temp_json, 'r') as f:
                            batch_data = json.load(f)
                            all_metadata.extend(batch_data if isinstance(batch_data, list) else [batch_data])
                        
                        os.remove(temp_json)
                    except subprocess.CalledProcessError as e:
                        logger.error(f"Batch {batch_num} failed: {e}")
                        if os.path.exists(temp_json): os.remove(temp_json)
                        for f in __import__('glob').glob(f"{config['getfq']['outdir']}/batch_*_temp.json"):
                            os.remove(f)
                        sys.exit(1)
                
                with open(json_path, 'w') as f:
                    json.dump(all_metadata, f, indent=2)
                logger.warning(f"Merged {num_batches} batches: {len(all_metadata)} entries")

        # Run getfq workflow
        profile = config['getfq']['profile']
        cmd = (
            " snakemake {profile} {dryrun} {jobs} {cores} --snakefile {snakefile} "
            " --directory {directory} --configfile {configfile} --nolock --latency-wait 300 "
            " --rerun-incomplete --rerun-triggers mtime --scheduler greedy "
        ).format(
            profile = f"--profile {profile}" if profile is not None else "",
            dryrun = f"--dry-run" if args.dryrun else "",
            jobs = f"--jobs {config['getfq']['jobs']}",
            cores = f"--cores {config['getfq']['jobs']}",
            snakefile = os.path.dirname(os.path.abspath(__file__))+"/workflow/rules/getfq.smk",
            directory = os.path.abspath(config['getfq']['outdir']),
            configfile = os.path.abspath(args.configfile)
        )

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            logger.error(e)
            sys.exit(1)


    if args.command == "runqc":
        
        # Check inputs and update configuration
        config = load_configfile(args.configfile)
        config = runqc_config_set(config, args)
        runqc_config_check(config)
        update_configfile(args.configfile, config)

        # Run runqc workflow
        profile = config['runqc']['profile']
        cmd = (
            " snakemake {profile} {dryrun} {jobs} {cores} --snakefile {snakefile} "
            " --directory {directory} --configfile {configfile} --nolock --latency-wait 60 "
            " --rerun-incomplete --rerun-triggers mtime --scheduler greedy "
        ).format(
            profile = f"--profile {profile}" if profile is not None else "",
            dryrun = f"--dry-run" if args.dryrun else "",
            jobs = f"--jobs {config['runqc']['jobs']}",
            cores = f"--cores {config['runqc']['jobs']}",
            snakefile = os.path.dirname(os.path.abspath(__file__))+"/workflow/rules/runqc.smk",
            directory = os.path.abspath(config['runqc']['outdir']),
            configfile = os.path.abspath(args.configfile)
        )

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            logger.error(e)
            sys.exit(1)
    

    if args.command == "run_rnaseq":
        
        # Check inputs and update configuration
        config = load_configfile(args.configfile)
        config = run_rnaseq_config_set(config, args)
        run_rnaseq_config_check(config)
        update_configfile(args.configfile, config)

        # Run run_rnaseq workflow
        profile = config['run_rnaseq']['profile']
        cmd = (
            " snakemake {profile} {dryrun} {jobs} {cores} --snakefile {snakefile} "
            " --directory {directory} --configfile {configfile} --nolock --latency-wait 60 "
            " --rerun-incomplete --rerun-triggers mtime --scheduler greedy "
        ).format(
            profile = f"--profile {profile}" if profile is not None else "",
            dryrun = f"--dry-run" if args.dryrun else "",
            jobs = f"--jobs {config['run_rnaseq']['jobs']}",
            cores = f"--cores {config['run_rnaseq']['jobs']}",
            snakefile = os.path.dirname(os.path.abspath(__file__))+"/workflow/rules/run_rnaseq.smk",
            directory = os.path.abspath(config['run_rnaseq']['outdir']),
            configfile = os.path.abspath(args.configfile)
        )

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            logger.error(e)
            sys.exit(1)


    if args.command == "run_atacseq":
        
        # Check inputs and update configuration
        config = load_configfile(args.configfile)
        config = run_atacseq_config_set(config, args)
        run_atacseq_config_check(config)
        update_configfile(args.configfile, config)

        # Run run_atacseq workflow
        profile = config['run_atacseq']['profile']
        cmd = (
            " snakemake {profile} {dryrun} {jobs} {cores} --snakefile {snakefile} "
            " --directory {directory} --configfile {configfile} --nolock --latency-wait 60 "
            " --rerun-incomplete --rerun-triggers mtime --scheduler greedy "
        ).format(
            profile = f"--profile {profile}" if profile is not None else "",
            dryrun = f"--dry-run" if args.dryrun else "",
            jobs = f"--jobs {config['run_atacseq']['jobs']}",
            cores = f"--cores {config['run_atacseq']['jobs']}",
            snakefile = os.path.dirname(os.path.abspath(__file__))+"/workflow/rules/run_atacseq.smk",
            directory = os.path.abspath(config['run_atacseq']['outdir']),
            configfile = os.path.abspath(args.configfile)
        )

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            logger.error(e)
            sys.exit(1)
