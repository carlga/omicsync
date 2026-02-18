import os
import sys
import re
import subprocess

from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d")

workflow_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
sys.path.append(os.path.dirname(workflow_dir))

from utils import load_json, find_key

# Set variables
GETREF_JSON = f"{config['getref']['species']}_ensembl_{config['getref']['release']}.json"
FILES = find_key(load_json(GETREF_JSON), 'ftp')
FILENAMES = [os.path.basename(file) for file in FILES]


rule all:
    input:
        expand(["CHECKSUMS_report.tsv", 
                "{filename}_extract_complete"], 
            filename=FILENAMES)


rule extract_ref_links:
    """
    Extracts download links of reference files.
    """
    input:
        json = GETREF_JSON
    output:
        tmp_links = temp(expand(["{filename}.tmp"], filename=FILENAMES))
    params:
        json = os.path.basename(GETREF_JSON)
    log: "logs/" + TIMESTAMP + "_extract_ref_links.log"
    threads: 1
    run:
        shell(" echo 'Extracting FTP links from '{params.json}'...' > {log} ")

        for file in FILES:
            name = os.path.basename(file)
            shell(" echo 'Found '{name}'!' ")
            shell(" echo {file} > {name}.tmp ")
        
        shell(" echo 'FTP links extracted!' >> {log} ")

rule download_ref:
    """
    Downloads reference files.
    """
    input:
        tmp_link = "{filename}.tmp"
    output:
        ref = "{filename}"
    log: "logs/" + TIMESTAMP + "_{filename}_download_ref.log"
    threads: 1
    run:
        with open(input.tmp_link, "r") as f:
            link = f.readline().strip()
        
        shell(" echo 'Downloading '{link}'...' > {log} ")

        shell(" curl {link} --remote-name --silent ")
        
        shell(" echo 'Finished downloading '{link}'' >> {log} ")

rule check_ref:
    """
    Evaluates reference file integrity.
    """
    input:
        ref = "{filename}"
    output:
        check = temp("{filename}_CHECKSUMS")
    params:
        ensembl_release = config['getref']['release']
    log: "logs/" + TIMESTAMP + "_{filename}_check_ref.log"
    threads: 1
    run:
        shell(" echo 'Checking '{input.ref}'...' > {log} ")

        for file in FILES:
            if input.ref in file:
                link = os.path.join(os.path.dirname(file), "CHECKSUMS")
                break

        checksum_origin = subprocess.check_output(
            f"curl {link} --silent", 
            shell=True, 
            text=True
        )
        if params.ensembl_release > 112:
            checksum_origin = re.search(
                r'^(\d+\s+\d+)\s+.*\s' + re.escape(input.ref) + r'$',
                checksum_origin,
                re.MULTILINE
            ).group(1)
        else:
            checksum_origin = re.search(
                r'^(\d+\s+\d+)\s+' + re.escape(input.ref) + r'$',
                checksum_origin,
                re.MULTILINE
            ).group(1)

        checksum_local = subprocess.check_output(
            f"sum {input.ref}", 
            shell=True, 
            text=True
        ).strip()

        test = "Passed" if checksum_origin == checksum_local else "Failed"
        result = f"{input.ref}\t{checksum_origin}\t{checksum_local}\t{test}"
        shell(" echo {result} > {output.check} ")
        shell(" echo 'Finished checking '{input.ref}'' >> {log} ")

rule check_ref_merge:
    """
    Combines reference file integrity reports.
    """
    input:
        checks = expand(["{filename}_CHECKSUMS"], filename=FILENAMES)
    output:
        check_report = "CHECKSUMS_report.tsv"
    log: "logs/" + TIMESTAMP + "_check_ref_merge.log"
    threads: 1
    shell:
        '''
        echo 'Generating file integrity report...' > {log}

        cat {input.checks} > {output.check_report}

        echo 'Finished report' >> {log}
        '''

rule extract_ref:
    """
    Uncompresses reference files.
    """
    input:
        ref = "{filename}",
        check_report = "CHECKSUMS_report.tsv"
    output:
        complete = temp("{filename}_extract_complete")
    log: "logs/" + TIMESTAMP + "_{filename}_extract_ref.log"
    threads: 1
    shell:
        '''
        echo 'Extracting '{input.ref}'...' > {log}

        gunzip {input.ref}

        echo 'Finished extracting '{input.ref}'' >> {log}

        touch {output.complete}
        '''


onsuccess:
    shell(" echo Download successful! ")
onerror:
    shell(" echo Download halted! ")
