"""
================================================================================
GETFQ MODULE - FASTQ File Download Workflow
================================================================================
Downloads FASTQ sequencing files from public repositories (ENA, SRA, etc.)
and validates file integrity using MD5 checksums.

Workflow steps:
1. Extract FTP download links from download metadata JSON
2. Download FASTQ files with robust retry logic
3. Validate file integrity using MD5 checksums
4. Generate integrity verification report

Supports:
- European Nucleotide Archive (ENA)
- NCBI Sequence Read Archive (SRA)
- Other FTP-accessible FASTQ repositories

Features:
- Automatic retry with exponential backoff for network failures
- Connection timeout protection (5 minutes per file)
- MD5 checksum validation against repository metadata
- Comprehensive error reporting
================================================================================
"""

import os
import sys
import re
import subprocess
from pathlib import Path

# Setup workflow environment
workflow_dir = Path(workflow.snakefile).parent.absolute()
sys.path.insert(0, str(workflow_dir.parent))

# Import utilities
from utils import load_json, find_key, filter_dict
from utils.workflow_utils import get_timestamp

# Module identifier
MODULE = 'getfq'
TIMESTAMP = get_timestamp()


# ================================================================================
# CONFIGURATION AND SETUP
# ================================================================================

# Load download metadata (silent)
GETFQ_JSON = "download_meta.json"
FILES = find_key(load_json(GETFQ_JSON), 'url')
FILES = [f for f in FILES if 'fastq.gz' in f and 'ftp' in f]
FILENAMES = [os.path.basename(f) for f in FILES]

# Validate that FASTQ files were found
if not FILENAMES:
    print(f"[{MODULE}] ✗ ERROR: No FASTQ files found in {GETFQ_JSON}!")
    print(f"[{MODULE}] Please check that the metadata file contains valid FTP URLs.")
    sys.exit(1)


# ================================================================================
# WILDCARD CONSTRAINTS
# ================================================================================

wildcard_constraints:
    filename=r"[^/]+\.fastq\.gz"


# ================================================================================
# WORKFLOW INITIALIZATION
# ================================================================================

onstart:
    # Detect sequencing layout
    paired_end = any('_1.fastq.gz' in f or '_2.fastq.gz' in f for f in FILENAMES)
    layout = "Paired-end" if paired_end else "Single-end"
    
    # Count samples (pairs count as one sample in PE)
    if paired_end:
        samples = len([f for f in FILENAMES if '_1.fastq.gz' in f])
    else:
        samples = len(FILENAMES)
    
    # Detect repository source from URLs
    if FILES:
        first_url = FILES[0]
        if 'ebi.ac.uk' in first_url or 'ena' in first_url:
            source = "ENA (European Nucleotide Archive)"
        elif 'ncbi' in first_url or 'sra' in first_url:
            source = "SRA (NCBI Sequence Read Archive)"
        else:
            source = "FTP repository"
    else:
        source = "Unknown"
    
    files_list = "\n".join([f"[{MODULE}]   • {f}" for f in FILENAMES])
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] FASTQ File Download Workflow'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Source: {source}'
        echo '[{MODULE}] Files: {len(FILENAMES)} ({samples} sample(s), {layout})'
        echo ''
        echo '[{MODULE}] Files to download:'
        echo '{files_list}'
        echo ''
        echo '[{MODULE}] Download settings:'
        echo '[{MODULE}]   • Retry attempts: 8 with 15s delays'
        echo '[{MODULE}]   • Connection timeout: 60s'
        echo '[{MODULE}]   • Max time per file: 60 min'
        echo '[{MODULE}]   • Speed limit threshold: 20 KB/s'
        echo '[{MODULE}]   • Resume support: enabled'
        echo ''
        echo '[{MODULE}] Workflow steps:'
        echo '[{MODULE}]   1. Extract FTP links from metadata'
        echo '[{MODULE}]   2. Download FASTQ files'
        echo '[{MODULE}]   3. Validate MD5 checksums'
        echo '[{MODULE}]   4. Generate verification report'
        echo ''
        echo '[{MODULE}] Timestamp: {TIMESTAMP}'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo ''
    """)


# ================================================================================
# TARGET RULES
# ================================================================================

rule all:
    input:
        expand("{filename}", filename=FILENAMES),
        "CHECKSUMS_report.tsv"


# ================================================================================
# FTP LINK EXTRACTION
# ================================================================================

rule extract_fq_links:
    """
    Extracts FTP download URLs from the metadata JSON file.
    Creates temporary files containing individual download links for each FASTQ file.
    Only processes files with 'fastq.gz' extension and 'ftp' protocol.
    """
    input:
        json = GETFQ_JSON
    output:
        tmp_links = temp(expand("{filename}.tmp", filename=FILENAMES))
    log: f"logs/{TIMESTAMP}_extract_fq_links.log"
    threads: 1
    run:
        shell(f"echo '[{MODULE}] Extracting FTP download links from metadata...' | tee {{log}}")
        
        for i, file in enumerate(FILES, 1):
            name = os.path.basename(file)
            shell(f"echo {file} > {name}.tmp")
            shell(f"echo '[{MODULE}]   [{i}/{len(FILES)}] Found {name}' | tee -a {{log}}")
        
        shell(f"echo '[{MODULE}] ✓ FTP links extracted for {len(FILES)} FASTQ files' | tee -a {{log}}")


# ================================================================================
# FASTQ FILE DOWNLOAD
# ================================================================================

rule download_fq:
    """
    Downloads FASTQ files from FTP servers with robust error handling.
    
    Network resilience features:
    - 8 retry attempts with 15-second delays between attempts
    - Automatic retry on connection refused errors
    - 60-minute timeout per file to prevent hanging
    - Resume support for interrupted downloads
    - Speed limit monitoring (fails if < 20 KB/s for 10 minutes)
    - Silent mode with error reporting only on failure
    
    Large FASTQ files may take considerable time depending on network speed.
    """
    input:
        tmp_link = "{filename}.tmp"
    output:
        fq = "{filename}"
    log: f"logs/{TIMESTAMP}_{{filename}}_download_fq.log"
    threads: 1
    run:
        with open(input.tmp_link, "r") as f:
            link = f.readline().strip()

        shell(f"""
            echo '[{MODULE}] Downloading {{wildcards.filename}}...' | tee {{log}}
            echo '[{MODULE}] Source: {link}' | tee -a {{log}}

            curl {link} \
                --remote-name \
                --silent \
                --show-error --fail-early \
                --retry 8 --retry-delay 15 --retry-connrefused --retry-max-time 1200 \
                --max-time 3600 --connect-timeout 60 \
                --speed-time 600 --speed-limit 20000 \
                --continue-at - \
                2>> {{log}} || {{{{
                echo '[{MODULE}] ✗ ERROR: Download failed for {{wildcards.filename}} after 8 attempts' | tee -a {{log}}
                echo '[{MODULE}] Check network connectivity and FTP server availability' | tee -a {{log}}
                exit 1
            }}}}

            file_size=$(du -h {{output.fq}} | cut -f1)
            echo '[{MODULE}] ✓ Successfully downloaded {{wildcards.filename}} ('$file_size')' | tee -a {{log}}
        """)


# ================================================================================
# INTEGRITY VERIFICATION
# ================================================================================

rule check_fq:
    """
    Validates FASTQ file integrity using MD5 checksum verification.
    Compares local MD5 hash against repository-provided checksum from metadata.
    
    MD5 checksums ensure that:
    - Files were downloaded completely without truncation
    - No corruption occurred during transfer
    - Files match exactly what the repository provides
    
    Generates per-file verification report with pass/fail status.
    """
    input:
        fq = "{filename}"
    output:
        check = temp("{filename}_CHECKSUMS")
    log: f"logs/{TIMESTAMP}_{{filename}}_check_fq.log"
    threads: 1
    run:
        shell(f"echo '[{MODULE}] Verifying integrity of {{input.fq}}...' | tee {{log}}")

        # Find file metadata from JSON
        file_meta = None
        for file in FILES:
            if input.fq in file:
                file_meta = filter_dict(load_json(GETFQ_JSON), ('url', file))
                break

        if not file_meta:
            shell(f"""
                echo '[{MODULE}] ✗ ERROR: Could not find metadata for {{input.fq}}' | tee -a {{log}}
                exit 1
            """)

        # Extract expected MD5 checksum from metadata
        checksum_origin = find_key(file_meta[0], 'md5')[0]
        
        shell(f"echo '[{MODULE}] Expected MD5: {checksum_origin}' | tee -a {{log}}")
        shell(f"echo '[{MODULE}] Computing local MD5 checksum...' | tee -a {{log}}")

        # Calculate local MD5 checksum
        checksum_local = subprocess.check_output(
            f"md5sum {input.fq}", shell=True, text=True
        ).strip().split()[0]

        shell(f"echo '[{MODULE}] Local MD5:    {checksum_local}' | tee -a {{log}}")

        # Compare checksums
        test = "Passed" if checksum_origin == checksum_local else "Failed"
        result = f"{input.fq}\t{checksum_origin}\t{checksum_local}\t{test}"

        shell(f"echo '{result}' > {{output.check}}")
        
        if test == "Passed":
            shell(f"echo '[{MODULE}] ✓ Integrity check PASSED for {{input.fq}}' | tee -a {{log}}")
        else:
            shell(f"echo '[{MODULE}] ✗ WARNING: Integrity check FAILED for {{input.fq}}' | tee -a {{log}}")
            shell(f"echo '[{MODULE}]   This file may be corrupted and should be re-downloaded' | tee -a {{log}}")


# ================================================================================
# REPORT GENERATION
# ================================================================================

rule check_fq_merge:
    """
    Combines individual MD5 verification results into a single TSV report.
    Provides summary statistics of passed and failed verifications.
    
    Report format:
    Filename | Original_MD5 | Local_MD5 | Status
    
    If any files fail verification, they should be re-downloaded before
    proceeding with downstream analysis to ensure data quality.
    """
    input:
        checks = expand("{filename}_CHECKSUMS", filename=FILENAMES)
    output:
        check_report = "CHECKSUMS_report.tsv"
    log: f"logs/{TIMESTAMP}_check_fq_merge.log"
    threads: 1
    shell:
        f'''
        echo '[{MODULE}] Generating integrity verification report...' | tee {{log}}
        
        echo -e "Filename\\tOriginal_MD5\\tLocal_MD5\\tStatus" > {{output.check_report}}
        cat {{input.checks}} >> {{output.check_report}}

        passed=$(grep -c "Passed" {{output.check_report}} || true)
        failed=$(grep -c "Failed" {{output.check_report}} || true)
        total=$((passed + failed))
        
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━' | tee -a {{log}}
        echo '[{MODULE}] MD5 INTEGRITY CHECK SUMMARY:' | tee -a {{log}}
        echo '[{MODULE}]   Total files: '$total | tee -a {{log}}
        echo '[{MODULE}]   ✓ Passed: '$passed' files' | tee -a {{log}}
        echo '[{MODULE}]   ✗ Failed: '$failed' files' | tee -a {{log}}
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━' | tee -a {{log}}
        
        if [ "$failed" -gt 0 ]; then
            echo '[{MODULE}] ⚠ WARNING: '$failed' file(s) failed MD5 verification!' | tee -a {{log}}
            echo '[{MODULE}] These files should be re-downloaded before analysis' | tee -a {{log}}
            echo '[{MODULE}] Check CHECKSUMS_report.tsv for details' | tee -a {{log}}
        else
            echo '[{MODULE}] ✓ All FASTQ files passed MD5 verification' | tee -a {{log}}
        fi
        '''


# ================================================================================
# WORKFLOW COMPLETION
# ================================================================================

onsuccess:
    # Detect layout and count samples
    paired_end = any('_1.fastq.gz' in f or '_2.fastq.gz' in f for f in FILENAMES)
    layout = "paired-end" if paired_end else "single-end"
    
    if paired_end:
        samples = len([f for f in FILENAMES if '_1.fastq.gz' in f])
    else:
        samples = len(FILENAMES)
    
    # Calculate total download size
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✓ FASTQ download complete!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Downloaded: {len(FILENAMES)} files ({samples} sample(s), {layout})'
        
        total_size=$(du -ch {' '.join(FILENAMES)} 2>/dev/null | grep total | cut -f1 || echo "N/A")
        echo '[{MODULE}] Total size: '$total_size
        
        echo ''
        echo '[{MODULE}] Verification report: CHECKSUMS_report.tsv'
        echo ''
        echo '[{MODULE}] Next steps:'
        echo '[{MODULE}]   - Review CHECKSUMS_report.tsv for integrity status'
        echo '[{MODULE}]   - Use runqc workflow for quality control'
        echo '[{MODULE}]   - Proceed with processing and analysis'
        echo ''
    """)

onerror:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✗ FASTQ download FAILED!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Common issues:'
        echo '[{MODULE}]   - Network connectivity problems'
        echo '[{MODULE}]   - FTP server temporarily unavailable or overloaded'
        echo '[{MODULE}]   - Invalid URLs in download_meta.json'
        echo '[{MODULE}]   - Insufficient disk space for large FASTQ files'
        echo '[{MODULE}]   - Download speed too slow (< 20 KB/s sustained)'
        echo '[{MODULE}]   - MD5 checksum mismatch (corrupted transfer)'
        echo ''
        echo '[{MODULE}] Troubleshooting steps:'
        echo '[{MODULE}]   1. Check logs/ directory for detailed error messages'
        echo '[{MODULE}]   2. Verify network connectivity and FTP server status'
        echo '[{MODULE}]   3. Check available disk space'
        echo '[{MODULE}]   4. Retry download - curl will resume incomplete files'
        echo '[{MODULE}]   5. Review CHECKSUMS_report.tsv if partially complete'
        echo '[{MODULE}]   6. Consider alternative data source if persistent failures'
        echo ''
        echo '[{MODULE}] Note: Downloads automatically retry up to 8 times'
        echo '[{MODULE}] Failed at timestamp: {TIMESTAMP}'
        echo ''
    """)
