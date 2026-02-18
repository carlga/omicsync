import os
import json
import re
import pandas as pd
import plotly.express as px


def load_json(file):
    """ Loads contents of json file to dictionary. """

    file = os.path.abspath(file)
  
    with open(file, "r") as f:
        d = json.load(f)
    
    return d

def find_key(d, key):
    """ Recursively searches for key in listed dictionary and returns values. """

    values = []

    if isinstance(d, list):
        for element in d:
            nested_values = find_key(element, key)
            for nested_value in nested_values:
                values.append(nested_value)
    else:
        for k,v in d.items():
            if k == key:
                values.append(v)
            elif isinstance(v, dict) or isinstance(v, list):
                nested_values = find_key(v, key)
                for nested_value in nested_values:
                    values.append(nested_value)
    
    return values

def filter_dict(d, key_value_pair):
    """ Recursively searches and filters listed dictionary for a given key-value pair. """

    d_filt = []

    if isinstance(d, list):
        for element in d:
            nested_dicts = filter_dict(element, key_value_pair)
            for nested_dict in nested_dicts:
                d_filt.append(nested_dict)
    else:
        for k,v in d.items():
            if k == key_value_pair[0] and v == key_value_pair[1]:
                d_filt.append(d)
            elif isinstance(v, dict) or isinstance(v, list):
                nested_dicts = filter_dict(v, key_value_pair)
                for nested_dict in nested_dicts:
                    d_filt.append(nested_dict)
    
    return d_filt

def check_strandedness(file, paired):
    """ Evaluates strandedness from RSeQC infer_experiment output. """

    with open(file, "r") as f:
        text = f.read()

    if paired:
        failed_fraction = re.search(r'failed to determine: (\d\.\d+)', text).group(1)
        fwd = re.search(r'\"1\+\+,1--,2\+-,2-\+\": (\d\.\d+)', text).group(1)
        rev = re.search(r'\"1\+-,1-\+,2\+\+,2--\": (\d\.\d+)', text).group(1)
    else:
        failed_fraction = re.search(r'failed to determine: (\d\.\d+)', text).group(1)
        fwd = re.search(r'\"\+\+,--\": (\d\.\d+)', text).group(1)
        rev = re.search(r'\+-,-\+\": (\d\.\d+)', text).group(1)
    
    ratio = float(fwd) / float(rev)

    if ratio > 10:
        strandedness = 'fwd'
    elif ratio < 0.1:
        strandedness = 'rev'
    else:
        strandedness = 'unstranded'

    return strandedness

def plot_biotype_composition(featurecounts_file, html_output_file):
    """ Generates interactive plot with biotype composition per sample. """

    df = pd.read_csv(featurecounts_file, sep='\t', skiprows=[0], index_col=[0])
    df = df.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)

    # Clean sample names
    samples = [re.sub(r'hisat2_aligned/([^\W_]+).bam', '\\1', x) for x in list(df.columns)]
    df.columns = samples

    # Calculate reads per million
    df = df/1e6

    df.index.name = 'Biotype'
    df = df.reset_index()

    df = df.melt(
        id_vars='Biotype', 
        var_name='Sample', 
        value_name='Reads'
    )

    # Biotype composition plot
    fig = px.bar(df, x="Reads", y="Sample", color="Biotype", template="seaborn")

    fig.update_layout(
        title="Biotype composition",
        xaxis_title="Read count (millions)",
        yaxis_title="",
        width=900,
        height=len(samples)*55+250,
        legend=dict(
            title_text="",
            font_size=8,
            orientation="h",
            yanchor="top",
            y=-1/len(samples),
            xanchor="left",
            x=0
        )
    )

    fig.write_html(html_output_file)
