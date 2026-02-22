"""Utility functions for osync workflows.

This module provides helper functions for workflow execution, path resolution,
command building, and argument parsing.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd
import plotly.express as px

from osync.logger import logger


def load_json(file: str) -> Dict[str, Any]:
    """Load contents of JSON file to dictionary.
    
    Args:
        file: Path to JSON file
        
    Returns:
        Dictionary containing JSON data
    """
    file_path = os.path.abspath(file)
    
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    return data


def find_key(data: Any, key: str) -> List[Any]:
    """Recursively search for key in nested dictionary/list and return values.
    
    Args:
        data: Dictionary or list to search
        key: Key name to find
        
    Returns:
        List of values found for the specified key
    """
    values = []
    
    if isinstance(data, list):
        for element in data:
            nested_values = find_key(element, key)
            values.extend(nested_values)
    elif isinstance(data, dict):
        for k, v in data.items():
            if k == key:
                values.append(v)
            elif isinstance(v, (dict, list)):
                nested_values = find_key(v, key)
                values.extend(nested_values)
    
    return values


def filter_dict(data: Any, key_value_pair: Tuple[str, Any]) -> List[Dict[str, Any]]:
    """Recursively search and filter nested dictionary for a given key-value pair.
    
    Args:
        data: Dictionary or list to filter
        key_value_pair: Tuple of (key, value) to match
        
    Returns:
        List of dictionaries matching the key-value pair
    """
    filtered = []
    
    if isinstance(data, list):
        for element in data:
            nested_dicts = filter_dict(element, key_value_pair)
            filtered.extend(nested_dicts)
    elif isinstance(data, dict):
        for k, v in data.items():
            if k == key_value_pair[0] and v == key_value_pair[1]:
                filtered.append(data)
            elif isinstance(v, (dict, list)):
                nested_dicts = filter_dict(v, key_value_pair)
                filtered.extend(nested_dicts)
    
    return filtered


def check_strandedness(file: str, paired: bool) -> str:
    """Evaluate strandedness from RSeQC infer_experiment output.
    
    Args:
        file: Path to RSeQC output file
        paired: Whether data is paired-end
        
    Returns:
        Strandedness type: 'fwd', 'rev', or 'unstranded'
    """
    with open(file, 'r') as f:
        text = f.read()
    
    if paired:
        fwd = re.search(r'"1\+\+,1--,2\+-,2-\+": (\d\.\d+)', text).group(1)
        rev = re.search(r'"1\+-,1-\+,2\+\+,2--": (\d\.\d+)', text).group(1)
    else:
        fwd = re.search(r'"\+\+,--": (\d\.\d+)', text).group(1)
        rev = re.search(r'\+-,-\+": (\d\.\d+)', text).group(1)
    
    ratio = float(fwd) / float(rev)
    
    if ratio > 10:
        return 'fwd'
    elif ratio < 0.1:
        return 'rev'
    else:
        return 'unstranded'


def plot_biotype_composition(featurecounts_file: str, html_output_file: str) -> None:
    """Generate interactive plot with biotype composition per sample.
    
    Args:
        featurecounts_file: Path to featureCounts output file
        html_output_file: Path for output HTML plot
    """
    df = pd.read_csv(featurecounts_file, sep='\t', skiprows=[0], index_col=[0])
    df = df.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
    
    # Clean sample names
    samples = [re.sub(r'hisat2_aligned/([^\W_]+).bam', r'\1', x) for x in df.columns]
    df.columns = samples
    
    # Calculate reads per million
    df = df / 1e6
    
    df.index.name = 'Biotype'
    df = df.reset_index()
    
    df = df.melt(
        id_vars='Biotype',
        var_name='Sample',
        value_name='Reads'
    )
    
    # Create biotype composition plot
    fig = px.bar(df, x='Reads', y='Sample', color='Biotype', template='seaborn')
    
    fig.update_layout(
        title='Biotype composition',
        xaxis_title='Read count (millions)',
        yaxis_title='',
        width=900,
        height=len(samples) * 55 + 250,
        legend=dict(
            title_text='',
            font_size=8,
            orientation='h',
            yanchor='top',
            y=-1 / len(samples),
            xanchor='left',
            x=0
        )
    )
    
    fig.write_html(html_output_file)


def get_workflow_dir() -> Path:
    """Get the path to the workflow rules directory.
    
    Returns:
        Path to the workflow/rules directory
    """
    # Go up 2 levels from utils.py: utils/ -> workflow/ -> osync/
    return Path(__file__).parent.parent / 'rules'


def build_snakemake_command(config: Dict[str, Any], module: str,
                            configfile: str, dryrun: bool = False) -> str:
    """Build Snakemake command for workflow execution.
    
    Args:
        config: Configuration dictionary
        module: Module name
        configfile: Path to configuration file
        dryrun: Whether to perform a dry run
        
    Returns:
        Complete Snakemake command string
    """
    from osync.main import MODULES
    
    module_config = config[module]
    profile = module_config.get('profile')
    jobs = module_config.get('jobs', 1)
    outdir = Path(module_config['outdir']).resolve()
    snakefile = get_workflow_dir() / MODULES[module]['snakefile']
    latency_wait = module_config.get('latency_wait', 60)
    
    cmd_parts = [
        'snakemake',
        f'--profile {profile}' if profile else '',
        '--dry-run' if dryrun else '',
        f'--jobs {jobs}',
        f'--cores {jobs}',
        f'--snakefile {snakefile}',
        f'--directory {outdir}',
        f'--configfile {Path(configfile).resolve()}',
        '--nolock',
        f'--latency-wait {latency_wait}',
        '--rerun-incomplete',
        '--rerun-triggers mtime',
        '--scheduler greedy'
    ]
    
    return ' '.join(part for part in cmd_parts if part)


def execute_workflow(module: str, config: Dict[str, Any],
                     configfile: str, dryrun: bool = False) -> None:
    """Execute a workflow for the specified module.
    
    Args:
        module: Module name
        config: Configuration dictionary
        configfile: Path to configuration file
        dryrun: Whether to perform a dry run
        
    Raises:
        SystemExit: If workflow execution fails
    """
    cmd = build_snakemake_command(config, module, configfile, dryrun)
    
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f'Workflow execution failed: {e}')
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('Workflow interrupted by user')
        sys.exit(130)


def handle_getref_workflow(config: Dict[str, Any], args: Any) -> None:
    """Handle getref workflow with its unique pre-processing requirements.
    
    Args:
        config: Configuration dictionary
        args: Parsed command-line arguments
    """
    # Initialize directory for download
    json_path = Path(config['getref']['outdir']) / (
        f"{config['getref']['species']}_ensembl_{config['getref']['release']}.json"
    )
    
    if json_path.exists():
        logger.warning(f"File exists '{json_path}'! Overwriting...")
    
    cmd = (
        f"gget ref {config['getref']['species']} "
        f"--which {config['getref']['which']} "
        f"--release {config['getref']['release']} "
        f"--out {json_path}"
    )
    
    subprocess.check_call(cmd, shell=True, stderr=subprocess.DEVNULL)
    
    # Execute workflow
    execute_workflow('getref', config, args.configfile, args.dryrun)


def handle_getfq_workflow(config: Dict[str, Any], args: Any) -> None:
    """Handle getfq workflow with its unique batch processing logic.
    
    Args:
        config: Configuration dictionary
        args: Parsed command-line arguments
    """
    outdir = Path(config['getfq']['outdir'])
    json_path = outdir / 'download_meta.json'
    batch_size = config['getfq'].get('batch_size', 1)
    
    if json_path.exists():
        logger.warning(f"File exists '{json_path}'! Using...")
    else:
        ids_list = config['getfq']['ids'].split()
        total_ids = len(ids_list)
        ftp_flag = '--ftp' if config['getfq'].get('ftp', False) else ''
        
        logger.warning(f'Fetching metadata for {total_ids} accession(s)')
        
        if batch_size == 1:
            # Single ffq call for all IDs
            cmd = f'ffq {config["getfq"]["ids"]} {ftp_flag} -o {json_path}'
            retry_cmd = f"""for i in {{1..5}}; do
                echo "Retrieving download links (attempt $i/5)"
                {cmd}
                sleep 5 && [[ -s {json_path} ]] && exit 0 || sleep 5
            done || (echo "ffq failed all retries" && exit 1)"""
            
            subprocess.check_call(retry_cmd, shell=True, stderr=subprocess.DEVNULL)
        else:
            # Batch processing
            all_metadata = []
            num_batches = (total_ids + batch_size - 1) // batch_size
            
            for i in range(0, total_ids, batch_size):
                batch = ids_list[i:i + batch_size]
                batch_num = i // batch_size + 1
                temp_json = outdir / f'batch_{batch_num}_temp.json'
                
                logger.warning(f'Batch {batch_num}/{num_batches}: {len(batch)} ID(s)')
                cmd = f"ffq {' '.join(batch)} {ftp_flag} -o {temp_json}"
                retry_cmd = f"""for i in {{1..5}}; do
                    echo "  Retrieving download links (attempt $i/5)"
                    {cmd}
                    sleep 5 && [[ -s {temp_json} ]] && exit 0 || sleep 5
                done || exit 1"""
                
                try:
                    subprocess.check_call(retry_cmd, shell=True, stderr=subprocess.DEVNULL)
                    
                    with open(temp_json, 'r') as f:
                        batch_data = json.load(f)
                        all_metadata.extend(
                            batch_data if isinstance(batch_data, list) else [batch_data]
                        )
                    
                    temp_json.unlink()
                except subprocess.CalledProcessError as e:
                    logger.error(f'Batch {batch_num} failed: {e}')
                    if temp_json.exists():
                        temp_json.unlink()
                    # Clean up all temp files
                    for temp_file in outdir.glob('batch_*_temp.json'):
                        temp_file.unlink()
                    sys.exit(1)
            
            with open(json_path, 'w') as f:
                json.dump(all_metadata, f, indent=2)
            logger.warning(f'Merged {num_batches} batches: {len(all_metadata)} entries')
    
    # Execute workflow
    execute_workflow('getfq', config, args.configfile, args.dryrun)
