import os
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
    if not output.exists():
        output.mkdir(parents = True)
     
    seqs_json = output / (tsv.stem + '_seqs.json')
    aiupred_npz = output / (tsv.stem + '_aiupred.npz')
    aiupred_p_npz = output / (tsv.stem + '_aiupred_p.npz')
    plddt_npz = output / (tsv.stem + '_plddt.npz') 
    plddt_vlp_npz = output / (tsv.stem + '_plddt_vlp.npz')

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

    # DataFrame
    gene, aiupred_col, plddt_col, chm_col, avg_col = 'Gene', f'AIUPred (> {aiupred_threshold})', f'pLDDT (< {plddt_threshold})', 'Counterharmonic mean', 'Average'
    df = pd.DataFrame({
        gene: list(gene_aiupred_p),
        aiupred_col: list(gene_aiupred_p.values()),
        plddt_col: list(gene_plddt_vlp.values()),
    })
    df[chm_col] = (df[aiupred_col] ** 2 + df[plddt_col] ** 2) / (df[aiupred_col] + df[plddt_col])
    df[avg_col] = (df[aiupred_col] + df[plddt_col])/ 2

    highlight = core.highlight(tsv)
    df_highlight = df[df['Gene'].isin(highlight)]
    
    # Bins
    labels = ['0–20%', '20–40%', '40–60%', '60–80%', '80–100%']
    bins_aiupred_disorder = pd.cut(df[aiupred_col], bins = [0, 20, 40, 60, 80, 100], labels = labels, include_lowest = True)
    bins_plddtp_disorder = pd.cut(df[plddt_col], bins = [0, 20, 40, 60, 80, 100], labels = labels, include_lowest = True)
    bins_chm_disorder = pd.cut(df[chm_col], bins = [0, 20, 40, 60, 80, 100], labels = labels, include_lowest = True)
    bins_avg_disorder = pd.cut(df[avg_col], bins = [0, 20, 40, 60, 80, 100], labels = labels, include_lowest = True)

    # Bar plots
    figure.bar_disorder_distribution(bins_aiupred_disorder, f'AIUPred (> {aiupred_threshold})')
    figure.bar_disorder_distribution(bins_plddtp_disorder, f'pLDDT (< {plddt_threshold})')
    figure.bar_disorder_distribution(bins_chm_disorder, 'Counterharmonic mean (pLDDT, AIUPred)')  
    figure.bar_disorder_distribution(bins_avg_disorder, 'Average (pLDDT, AIUPred)')  

    # Scatter
    figure.scatter_disorder_distribution(df, df_highlight, bins_chm_disorder)
    figure.scatter_disorder_distribution(df, df_highlight, bins_avg_disorder)

if __name__ == '__main__':
    app()