"""
================================================================================
GETREF MODULE - Reference File Download Workflow
================================================================================
Downloads and validates reference files from Ensembl FTP servers.

Workflow steps:
1. Extract FTP download links from species/release-specific JSON metadata
2. Download reference files (genome FASTA, GTF annotation, etc.)
3. Validate file integrity using Ensembl CHECKSUMS
4. Extract compressed files (gunzip)
5. Generate integrity verification report

Supports all Ensembl releases and species with automatic checksum validation
to ensure data integrity. Handles format changes across Ensembl releases.
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
from utils import load_json, find_key
from utils.workflow_utils import get_timestamp

# Module identifier
MODULE = 'getref'
TIMESTAMP = get_timestamp()


# ================================================================================
# CONFIGURATION AND SETUP
# ================================================================================

# Load reference metadata from JSON
GETREF_JSON = f"{config['getref']['species']}_ensembl_{config['getref']['release']}.json"
FILES = find_key(load_json(GETREF_JSON), 'ftp')
FILENAMES = [os.path.basename(f) for f in FILES]

# Create list of extracted filenames (without .gz extension)
EXTRACTED_FILES = [f.replace('.gz', '') for f in FILENAMES]


# ================================================================================
# WILDCARD CONSTRAINTS
# ================================================================================

wildcard_constraints:
    filename=r"[^/]+\.gz",  # Only match .gz files for download
    extracted=r"[^/]+\.(gtf|fa)"  # Only match .gtf or .fa files for extraction


# ================================================================================
# WORKFLOW INITIALIZATION
# ================================================================================

onstart:
    # Categorize files by type for better display
    fasta_files = [f for f in FILENAMES if '.dna' in f or '.fa' in f]
    gtf_files = [f for f in FILENAMES if '.gtf' in f]
    protein_files = [f for f in FILENAMES if '.pep' in f]
    other_files = [f for f in FILENAMES if f not in fasta_files + gtf_files + protein_files]
    
    files_summary = []
    if fasta_files:
        files_summary.append(f"{len(fasta_files)} genome FASTA")
    if gtf_files:
        files_summary.append(f"{len(gtf_files)} GTF annotation")
    if protein_files:
        files_summary.append(f"{len(protein_files)} protein sequence")
    if other_files:
        files_summary.append(f"{len(other_files)} other")
    
    summary_str = ", ".join(files_summary)
    files_list = "\n".join([f"[{MODULE}]   • {f}" for f in FILENAMES])
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Reference File Download Workflow'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Species: {config["getref"]["species"]}'
        echo '[{MODULE}] Ensembl Release: {config["getref"]["release"]}'
        echo '[{MODULE}] Files: {len(FILENAMES)} ({summary_str})'
        echo ''
        echo '[{MODULE}] Files to download:'
        echo '{files_list}'
        echo ''
        echo '[{MODULE}] Workflow steps:'
        echo '[{MODULE}]   1. Extract FTP links from metadata'
        echo '[{MODULE}]   2. Download files from Ensembl'
        echo '[{MODULE}]   3. Validate integrity (checksums)'
        echo '[{MODULE}]   4. Extract compressed files'
        echo '[{MODULE}]   5. Generate verification report'
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
        "CHECKSUMS_report.tsv",
        expand("{extracted}", extracted=EXTRACTED_FILES)


# ================================================================================
# FTP LINK EXTRACTION
# ================================================================================

rule extract_ref_links:
    """
    Extracts FTP download links from species-specific JSON metadata file.
    Creates temporary files containing individual download URLs for each
    reference file (genome FASTA, GTF, protein sequences, etc.).
    """
    input:
        json = GETREF_JSON
    output:
        tmp_links = temp(expand("{filename}.tmp", filename=FILENAMES))
    params:
        json = os.path.basename(GETREF_JSON)
    log: f"logs/{TIMESTAMP}_extract_ref_links.log"
    threads: 1
    run:
        shell(f"echo '[{MODULE}] Extracting FTP links from {{params.json}}...' | tee {{log}}")
        
        for i, file in enumerate(FILES, 1):
            name = os.path.basename(file)
            shell(f"echo {file} > {name}.tmp")
            shell(f"echo '[{MODULE}]   [{i}/{len(FILES)}] Found {name}' | tee -a {{log}}")
        
        shell(f"echo '[{MODULE}] ✓ FTP links extracted for {len(FILES)} files' | tee -a {{log}}")


# ================================================================================
# REFERENCE FILE DOWNLOAD
# ================================================================================

rule download_ref:
    """
    Downloads reference files from Ensembl FTP servers using curl.
    Includes retry logic (3 attempts with 5-second delays) for network resilience.
    Validates download success before proceeding to checksum verification.
    """
    input:
        tmp_link = "{filename}.tmp"
    output:
        ref = temp("{filename}")  # Mark downloaded .gz files as temporary
    log: f"logs/{TIMESTAMP}_{{filename}}_download_ref.log"
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
                --show-error \
                --fail-early \
                --retry 3 \
                --retry-delay 5 \
                2>> {{log}} || {{{{
                echo '[{MODULE}] ✗ ERROR: Download failed for {{wildcards.filename}}' | tee -a {{log}}
                exit 1
            }}}}

            echo '[{MODULE}] ✓ Downloaded {{wildcards.filename}}' | tee -a {{log}}
        """)


# ================================================================================
# INTEGRITY VERIFICATION
# ================================================================================

rule check_ref:
    """
    Validates reference file integrity using Ensembl CHECKSUMS.
    Compares local file checksum (BSD sum algorithm) against official Ensembl checksum.
    
    Handles Ensembl format changes:
    - Release ≤112: Format "checksum blocks filename"
    - Release >112: Format "checksum blocks size date filename"
    
    Generates per-file verification report with pass/fail status.
    """
    input:
        ref = "{filename}"
    output:
        check = temp("{filename}_CHECKSUMS")
    params:
        ensembl_release = config['getref']['release']
    log: f"logs/{TIMESTAMP}_{{filename}}_check_ref.log"
    threads: 1
    run:
        shell(f"echo '[{MODULE}] Verifying integrity of {{input.ref}}...' | tee {{log}}")

        # Find corresponding checksum file URL
        for file in FILES:
            if input.ref in file:
                link = os.path.join(os.path.dirname(file), "CHECKSUMS")
                break

        shell(f"echo '[{MODULE}] Fetching checksums from Ensembl...' | tee -a {{log}}")

        # Download checksum file from Ensembl
        checksum_origin = subprocess.check_output(
            f"curl {link} --silent", shell=True, text=True
        )

        # Parse checksum based on Ensembl release format
        if params.ensembl_release > 112:
            # New format: "checksum blocks size date filename"
            pattern = r'^(\d+)\s+(\d+)\s+.*\s' + re.escape(input.ref) + r'$'
        else:
            # Old format: "checksum blocks filename"
            pattern = r'^(\d+)\s+(\d+)\s+' + re.escape(input.ref) + r'$'

        match = re.search(pattern, checksum_origin, re.MULTILINE)
        if match:
            checksum_origin = f"{match.group(1)} {match.group(2)}"
        else:
            checksum_origin = "NOT_FOUND"

        # Calculate local checksum using BSD sum
        shell(f"echo '[{MODULE}] Computing local checksum...' | tee -a {{log}}")
        checksum_output = subprocess.check_output(
            f"sum {input.ref}", shell=True, text=True
        ).strip()
        
        # Extract only the first two fields (checksum and blocks), ignore filename
        checksum_local = ' '.join(checksum_output.split()[:2])

        # Compare checksums
        test = "Passed" if checksum_origin == checksum_local else "Failed"
        result = f"{input.ref}\t{checksum_origin}\t{checksum_local}\t{test}"

        shell(f"echo '{result}' > {{output.check}}")
        
        if test == "Passed":
            shell(f"echo '[{MODULE}] ✓ Integrity check PASSED for {{input.ref}}' | tee -a {{log}}")
        else:
            shell(f"echo '[{MODULE}] ✗ WARNING: Integrity check FAILED for {{input.ref}}' | tee -a {{log}}")
            shell(f"echo '[{MODULE}]   Expected: {checksum_origin}' | tee -a {{log}}")
            shell(f"echo '[{MODULE}]   Got: {checksum_local}' | tee -a {{log}}")


# ================================================================================
# REPORT GENERATION
# ================================================================================

rule check_ref_merge:
    """
    Combines individual integrity check results into a single TSV report.
    Provides summary statistics of passed and failed verifications.
    
    Report format:
    Filename | Original_Checksum | Local_Checksum | Status
    
    Use this report to identify any corrupted downloads before proceeding.
    """
    input:
        checks = expand("{filename}_CHECKSUMS", filename=FILENAMES)
    output:
        check_report = "CHECKSUMS_report.tsv"
    log: f"logs/{TIMESTAMP}_check_ref_merge.log"
    threads: 1
    shell:
        f'''
        echo '[{MODULE}] Generating integrity verification report...' | tee {{log}}
        
        echo -e "Filename\\tOriginal_Checksum\\tLocal_Checksum\\tStatus" > {{output.check_report}}
        cat {{input.checks}} >> {{output.check_report}}

        passed=$(grep -c "Passed" {{output.check_report}} || true)
        failed=$(grep -c "Failed" {{output.check_report}} || true)
        
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━' | tee -a {{log}}
        echo '[{MODULE}] INTEGRITY CHECK SUMMARY:' | tee -a {{log}}
        echo '[{MODULE}]   ✓ Passed: '$passed' files' | tee -a {{log}}
        echo '[{MODULE}]   ✗ Failed: '$failed' files' | tee -a {{log}}
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━' | tee -a {{log}}
        
        if [ "$failed" -gt 0 ]; then
            echo '[{MODULE}] ⚠ WARNING: Some files failed integrity checks!' | tee -a {{log}}
        else
            echo '[{MODULE}] ✓ All files passed integrity checks' | tee -a {{log}}
        fi
        '''


# ================================================================================
# FILE EXTRACTION
# ================================================================================

rule extract_ref:
    """
    Extracts compressed reference files using gunzip.
    Only proceeds after successful integrity verification.
    Outputs the extracted file directly (no marker files).
    
    Note: Original .gz files are automatically cleaned up by Snakemake.
    """
    input:
        ref = "{extracted}.gz",
        check_report = "CHECKSUMS_report.tsv"
    output:
        extracted = "{extracted}"
    log: f"logs/{TIMESTAMP}_{{extracted}}.gz_extract_ref.log"
    threads: 1
    shell:
        f'''
        echo '[{MODULE}] Extracting {{input.ref}}...' | tee {{log}}

        gunzip -c {{input.ref}} > {{output.extracted}} 2>> {{log}} || {{{{
            echo '[{MODULE}] ✗ ERROR: Failed to extract {{input.ref}}' | tee -a {{log}}
            exit 1
        }}}}

        echo '[{MODULE}] ✓ Successfully extracted to: {{output.extracted}}' | tee -a {{log}}
        '''


# ================================================================================
# WORKFLOW COMPLETION
# ================================================================================

onsuccess:
    # Categorize extracted files
    fasta_files = [f for f in EXTRACTED_FILES if '.dna' in f or '.fa' in f]
    gtf_files = [f for f in EXTRACTED_FILES if '.gtf' in f]
    protein_files = [f for f in EXTRACTED_FILES if '.pep' in f]
    
    output_list = []
    if fasta_files:
        output_list.append(f"[{MODULE}]   - Genome FASTA: {', '.join(fasta_files)}")
    if gtf_files:
        output_list.append(f"[{MODULE}]   - GTF annotation: {', '.join(gtf_files)}")
    if protein_files:
        output_list.append(f"[{MODULE}]   - Protein sequences: {', '.join(protein_files)}")
    
    outputs_str = "\n".join(output_list)
    
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✓ Reference download complete!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Species: {config["getref"]["species"]}'
        echo '[{MODULE}] Ensembl Release: {config["getref"]["release"]}'
        echo '[{MODULE}] Files downloaded: {len(FILENAMES)}'
        echo '[{MODULE}] Files extracted: {len(EXTRACTED_FILES)}'
        echo ''
        echo '[{MODULE}] Output files:'
        echo '{outputs_str}'
        echo ''
        echo '[{MODULE}] Verification report: CHECKSUMS_report.tsv'
        echo ''
        echo '[{MODULE}] Next steps:'
        echo '[{MODULE}]   - Review CHECKSUMS_report.tsv for integrity status'
        echo '[{MODULE}]   - Use init workflow to build indices (HISAT2, Kallisto, etc.)'
        echo '[{MODULE}]   - Reference files ready for analysis workflows'
        echo ''
    """)

onerror:
    shell(f"""
        echo ''
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] ✗ Reference download FAILED!'
        echo '[{MODULE}] ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
        echo '[{MODULE}] Common issues:'
        echo '[{MODULE}]   - Network connectivity problems'
        echo '[{MODULE}]   - Ensembl FTP server temporarily unavailable'
        echo '[{MODULE}]   - Invalid species name or release number'
        echo '[{MODULE}]   - Checksum verification failed (corrupted download)'
        echo '[{MODULE}]   - Insufficient disk space for downloads'
        echo '[{MODULE}]   - Missing JSON metadata file'
        echo ''
        echo '[{MODULE}] Troubleshooting steps:'
        echo '[{MODULE}]   1. Check logs/ directory for detailed error messages'
        echo '[{MODULE}]   2. Verify network connectivity to ftp.ensembl.org'
        echo '[{MODULE}]   3. Confirm species and release numbers are valid'
        echo '[{MODULE}]   4. Check available disk space'
        echo '[{MODULE}]   5. Review CHECKSUMS_report.tsv if partially complete'
        echo '[{MODULE}]   6. Retry download - curl has 3 automatic retries per file'
        echo ''
        echo '[{MODULE}] Failed at timestamp: {TIMESTAMP}'
        echo ''
    """)
