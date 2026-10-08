from pathlib import Path
from typing import Annotated
import json

import typer
import numpy as np
import pandas as pd

from . import core
from . import figure


app = typer.Typer()


@app.command()
def compare(tsv: Annotated[Path, typer.Argument(exists = True, dir_okay = False, help = 'Path to the TSV with Gene and Uniprot columns.')],
            output: Annotated[Path, typer.Argument(dir_okay = True, help = 'Path to output folder')],
            aiupred_threshold: Annotated[float, typer.Option('-a', '--aiupred-threshold', help = 'AIUPred threshold (disorder if above).')] = 0.5,
            plddt_threshold: Annotated[int, typer.Option('-p', '--plddt-threshold', help = 'pLDDT threshold (disorder if below).')] = 80):
    output_cache = output / 'cache'
    output_equal_bins = output / 'equal_bins'
    output_ordered_moderate_disordered_bins = output / 'ordered_moderate_disordered_bins'    
    if not output.exists():
        output.mkdir(parents = True)
        output_cache.mkdir()
        output_equal_bins.mkdir()
        output_ordered_moderate_disordered_bins.mkdir()
     
    seqs_json = output_cache / (tsv.stem + '_seqs.json')
    aiupred_npz = output_cache / (tsv.stem + '_aiupred.npz')
    aiupred_p_npz = output_cache / (tsv.stem + '_aiupred_p.npz')
    plddt_npz = output_cache / (tsv.stem + '_plddt.npz') 
    plddt_vlp_npz = output_cache / (tsv.stem + '_plddt_vlp.npz')

    '''
    Profiler
    '''
    typer.echo('Retrieving UniProt identifiers...')
    gene_uid = core.gene_to_uid(tsv)
    typer.echo(f'└── {len(gene_uid)} UniProt identifiers were obtained.\n')

    if not seqs_json.exists():
        typer.echo('Downloading protein sequences from the UniProt database...')
        gene_seq = core.gene_to_seq(gene_uid)
        with open(seqs_json, 'w') as output:
            json.dump(gene_seq, output)
        typer.echo(f'└── Downloaded {len(gene_seq)} protein sequences from UniProt \n')
    else:
        typer.echo('Cache detected for protein sequences. Loading from cache...')
        with open(seqs_json, 'rt') as input:
            gene_seq = json.load(input)
        typer.echo(f'└── Loaded {len(gene_seq)} protein sequences from cache \n')

    if not aiupred_npz.exists():
        typer.echo('Processing AIUPred...')
        gene_aiupred = core.gene_to_aiupred(gene_seq)
        np.savez(aiupred_npz, **gene_aiupred)
        typer.echo(f'└── Processed {len(gene_aiupred)} protein sequences\n')
    else:
        typer.echo('AIUPred cache detected. Loading from cache...')
        loaded = np.load(aiupred_npz)
        gene_aiupred = {gene: loaded[gene] for gene in loaded}
        typer.echo(f'└── Loaded AIUPred for {len(gene_aiupred)} protein sequences from cache\n')

    if not aiupred_p_npz.exists():
        typer.echo('Calculating disorder percentage based on AIUPred...')
        gene_aiupred_p = core.gene_to_aiupred_p(gene_aiupred, aiupred_threshold)
        np.savez(aiupred_p_npz, **gene_aiupred_p)
        typer.echo(f'└── Calculated disorder percentage based on AIUPred for {len(gene_aiupred_p)} proteins\n')
    else:
        typer.echo('Disorder percentage based on AIUPred cache detected. Loading from cache...')
        loaded = np.load(aiupred_p_npz)
        gene_aiupred_p = {gene: loaded[gene] for gene in loaded}
        typer.echo(f'└── Loaded disorder percentage based on AIUPred for {len(gene_aiupred_p)} proteins from cache\n')

    if not plddt_npz.exists():
        typer.echo('Downloading pLDDT from AlphaFold...')
        gene_plddt = core.gene_to_plddt(gene_uid)
        np.savez(plddt_npz, **gene_plddt)
        typer.echo(f'└── Downloaded pLDDT for {len(gene_plddt)} proteins from AlphaFold\n')
    else:
        typer.echo('pLDDT cache detected. Loading from cache...')
        loaded = np.load(plddt_npz)
        gene_plddt = {gene: loaded[gene] for gene in loaded}
        typer.echo(f'└── Loaded per-residue pLDDT for {len(gene_plddt)} proteins from cache\n')

    if not plddt_vlp_npz.exists():
        typer.echo('Calculating disorder percentage based on pLDDT...')
        gene_plddt_vlp = core.gene_to_plddt_vlp(gene_plddt, plddt_threshold)
        np.savez(plddt_vlp_npz, **gene_plddt_vlp)
        typer.echo(f'└── Calculated disorder percentage based on pLDDT for {len(gene_plddt_vlp)} proteins\n')
    else:
        typer.echo('Disorder percentage based on pLDDT cache detected. Loading from cache...')
        loaded = np.load(plddt_vlp_npz)
        gene_plddt_vlp = {gene: loaded[gene] for gene in loaded}
        typer.echo(f'└── Loaded disorder percentage based on pLDDT for {len(gene_plddt_vlp)} proteins from cache\n')
    
    '''
    Figure
    '''
    common = set(gene_aiupred_p) & set(gene_plddt_vlp)
    
    gene_aiupred_p = {g : gene_aiupred_p[g] * 100  for g in common}
    gene_plddt_vlp = {g : gene_plddt_vlp[g] * 100 for g in common}

    # Data
        # Columns
    GENE_COL = 'Gene'
    AIUPRED_COL = f'PPIDR-AIUPred (residues with score < {aiupred_threshold} were considered disordered)'
    PLDDT_COL = f'PPIDR-pLDDT (residues with score < {plddt_threshold} were considered disordered)'
    CHM_COL = 'Contraharmonic mean between PPIDR-AIUpred and PPIDR-pLDDT'
    AVG_COL = 'Average between PPIDR-AIUpred and PPIDR-pLDDT'

    PREDICTOR_COLS = (AIUPRED_COL, PLDDT_COL)
    AGGREGATE_COLS = (AVG_COL, CHM_COL)
    PPIDR_COLS = PREDICTOR_COLS + AGGREGATE_COLS
    ALL_COLLS = (GENE_COL, ) + PPIDR_COLS

        # DataFrame 
    df = pd.DataFrame()
    df[GENE_COL] = list(gene_aiupred_p)
    df[AIUPRED_COL] = list(gene_aiupred_p.values())
    df[PLDDT_COL] = list(gene_plddt_vlp.values())
    df[CHM_COL] = (df[AIUPRED_COL] ** 2 + df[PLDDT_COL] ** 2) / (df[AIUPRED_COL] + df[PLDDT_COL])
    df[AVG_COL] = (df[AIUPRED_COL] + df[PLDDT_COL])/ 2

    df['PR_AIUPRED'] = df[AIUPRED_COL].rank(pct = True, method = 'max') * 100
    df['PR_PLDDT_COL'] = df[PLDDT_COL].rank(pct = True, method = 'max') * 100
    df['PR_AVG_COL'] = df[AVG_COL].rank(pct = True, method = 'max') * 100
    df['PR_CHM_COL'] = df[CHM_COL].rank(pct = True, method = 'max') * 100
    
        # Highlight DataFrame
    highlight = core.highlight(tsv)
    targets = df[df['Gene'].isin(highlight)]

        # Save
    df.to_csv(output / 'summary.csv', index = False)
    targets.to_csv(output / 'targets_summary.csv', index = False)

    # Correlation PPIDR-AIUPred and PPIDR-pLDDT
    corrs = core.corr(list(gene_aiupred_p.values()), list(gene_plddt_vlp.values()))
    corrs_df = pd.DataFrame(corrs)
    corrs_df.to_csv(output / 'corrs.csv', index = False)

    # Binning
    BINNING_SCHEMA = {
        'EQUAL': {
            'bins' : [0, 20, 40, 60, 80, 100],
            'labels': ['0-20% disorder', '20-40% disorder', '40-60% disorder', '60-80% disorder', '80-100% disorder'],
            'bar_colors': ["#2B5375", "#A8BEDC", "#C6B4D8", "#E8B4C0", "#E18181"],
            'scatter_colors': ["#2B5375", "#A8BEDC", "#409151", "#E8B4C0", "#EB3B3B"],
            'output': output_equal_bins
        },
        'ORDERED_MODERATE_DISORDERED': {
            'bins': [0, 10, 30, 100],
            'labels': ['Structured proteins (< 10% disorder)', 'Moderately disordered proteins (20–40% disorder)', 'Highly disordered proteins (> 30% disorder)'],
            'bar_colors': ["#2B5375", "#A8BEDC", "#E18181"],
            'scatter_colors': ["#2B5375", "#789AD1", "#EB3B3B"],
            'output': output_ordered_moderate_disordered_bins
        }
    }

    for schema in BINNING_SCHEMA.values():
        # Bar plot
        for col in PPIDR_COLS:
            bins = pd.cut(df[col], bins = schema['bins'], labels = schema['labels'], right = False, include_lowest = True)
            fig = figure.bar_disorder_distribution(bins, schema['bar_colors'],col)
            fig.write_html(schema['output'] / f'BAR_{col}.html')

        # Scatter plot
        for col in AGGREGATE_COLS:
            bins = pd.cut(df[col], bins = schema['bins'], labels = schema['labels'], include_lowest = True)
            bins = bins.cat.set_categories(schema['labels'], ordered=True)
            fig = figure.scatter_disorder_distribution(df, AIUPRED_COL, PLDDT_COL, bins, schema['scatter_colors'], col, targets, ALL_COLLS)
            fig.write_html(schema['output'] / f'SCATTER_{col}.html')

        # Interceptions
        compare_bins = []
        for col in AGGREGATE_COLS:
            bins = pd.cut(df[col], bins = schema['bins'], labels = schema['labels'], include_lowest = True)
            compare_bins.append(bins)
        fig = figure.sankey_distribution_compare(compare_bins, AGGREGATE_COLS)
        fig.write_html(schema['output'] / f'SANKEY.html')

if __name__ == '__main__':
    app()