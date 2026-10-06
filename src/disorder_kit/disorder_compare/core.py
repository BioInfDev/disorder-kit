from pathlib import Path

import typer
import numpy as np
from scipy.stats import pearsonr, spearmanr
from pandas import read_csv
from aiupred import AIUPred

from .. import fetcher


def load_npz_to_dict(path: Path) -> dict:
    with np.load(path) as data:
        return {i: data[i] for i in data.files}

def gene_to_uid(tsv: Path) -> dict[str, str]:
    tsv = read_csv(tsv, sep = '\t').dropna(subset = ['Uniprot'])
    return {gene: uid for gene, uid in zip(tsv['Gene'], tsv['Uniprot'])}

def highlight(tsv: Path):
    tsv = read_csv(tsv, sep = '\t').dropna(subset = ['Gene'])
    return [gene for gene, highlight in zip(tsv['Gene'], tsv['Highlight']) if highlight == '+']

def gene_to_seq(gene_uid: dict) -> dict[str, list]:
    with typer.progressbar(gene_uid.items(), length = len(gene_uid)) as progress:
        return {gene: fetcher.get_protein_sequence(uid) for gene, uid in progress}

def gene_to_aiupred(gene_seq: dict) -> dict[str, np.ndarray]:
    predictor = AIUPred()
    with typer.progressbar(gene_seq.items(), length = len(gene_seq)) as progress:
        return {gene: predictor.predict_disorder(seq) for gene, seq in progress}

def gene_to_aiupred_p(gene_aiupred: dict, threshold) -> dict[str, float]:
    with typer.progressbar(gene_aiupred.items(), length = len(gene_aiupred)) as progress:
        return {gene: sum(x > threshold for x in aiupred) / len(aiupred) for gene, aiupred in progress}

def gene_to_plddt(gene_uid: dict) -> dict[str, np.ndarray]:
    with typer.progressbar(gene_uid.items(), length = len(gene_uid)) as progress:
        return {gene: plddt for gene, uid in progress if (plddt := fetcher.get_pLDDT(uid)) is not None and len(plddt) > 0}

def gene_to_plddt_vlp(gene_plddt: dict, threshold) -> dict[str, float]:
    with typer.progressbar(gene_plddt.items(), length = len(gene_plddt)) as progress:
        return {gene: sum(x < threshold for x in plddt) / len(plddt) for gene, plddt in progress}

def corr(l1: list, l2: list) -> dict[str, list]:
    r_p, p_p = pearsonr(l1, l2)
    r_s, p_s = spearmanr(l1, l2)
    return {'Pearson': [r_p, p_p], 'Spearman': [r_s, p_s]}
