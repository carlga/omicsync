import os
import sys
import re
import subprocess

from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d")

workflow_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
sys.path.append(os.path.dirname(workflow_dir))

from utils import load_json, find_key, filter_dict

# Set variables
GETFQ_JSON = "download_meta.json"
FILES = find_key(load_json(GETFQ_JSON), 'url')
FILES = [file for file in FILES if 'fastq.gz' in file and 'ftp' in file]
FILENAMES = [os.path.basename(file) for file in FILES]

if not FILENAMES:
    sys.exit("No fastq files found!")


rule all:
    input:
        expand(["{filename}",
                "CHECKSUMS_report.tsv"], 
            filename=FILENAMES)


rule extract_fq_links:
    """
    Extracts download links of fastq files.
    """
    input:
        json = GETFQ_JSON
    output:
        tmp_links = temp(expand(["{filename}.tmp"], filename=FILENAMES))
    log: "logs/" + TIMESTAMP + "_extract_fq_links.log"
    threads: 1
    run:
        shell(" echo 'Extracting FTP links...' > {log} ")

        for file in FILES:
            name = os.path.basename(file)
            shell(" echo {file} > {name}.tmp ")
        
        shell(" echo 'FTP links extracted!' >> {log} ")

rule download_fq:
    """
    Downloads fastq files.
    """
    input:
        tmp_link = "{filename}.tmp"
    output:
        fq = "{filename}"
    log: "logs/" + TIMESTAMP + "_{filename}_download_fq.log"
    threads: 1
    run:
        with open(input.tmp_link, "r") as f:
            link = f.readline().strip()
        
        shell("""
            echo 'Downloading {link}...' | tee -a {log}
            
            curl {link} \
                --remote-name \
                --silent \
                --fail-early \
                --retry 5 \
                --retry-delay 10 \
                --retry-connrefused \
                --max-time 300 || {{
                echo 'FAILED after 5 retries: {link}' | tee -a {log}
                exit 1
            }}
            
            echo 'Successfully downloaded {link}' | tee -a {log}
        """)

rule check_fq:
    """
    Evaluates fastq file integrity.
    """
    input:
        fq = "{filename}"
    output:
        check = temp("{filename}_CHECKSUMS")
    log: "logs/" + TIMESTAMP + "_{filename}_check_fq.log"
    threads: 1
    run:
        shell(" echo 'Checking '{input.fq}'...' > {log} ")

        for file in FILES:
            if input.fq in file:
                file_meta = filter_dict(load_json(GETFQ_JSON), ('url', file))
                break
        
        checksum_origin = find_key(file_meta[0], 'md5')[0]
        
        checksum_local = subprocess.check_output(
            f"md5sum {input.fq}", 
            shell=True, 
            text=True
        ).strip()
        checksum_local = re.split(r'\s+', checksum_local)[0]

        test = "Passed" if checksum_origin == checksum_local else "Failed"
        result = f"{input.fq}\t{checksum_origin}\t{checksum_local}\t{test}"
        shell(" echo {result} > {output.check} ")
        shell(" echo 'Finished checking '{input.fq}'' >> {log} ")

rule check_fq_merge:
    """
    Combines fastq file integrity reports.
    """
    input:
        checks = expand(["{filename}_CHECKSUMS"], filename=FILENAMES)
    output:
        check_report = "CHECKSUMS_report.tsv"
    log: "logs/" + TIMESTAMP + "_check_fq_merge.log"
    threads: 1
    shell:
        '''
        echo 'Generating fq file integrity report...' > {log}

        cat {input.checks} > {output.check_report}

        echo 'Finished report' >> {log}
        '''


onsuccess:
    shell(" echo Download successful! ")
onerror:
    shell(" echo Download halted! ")
