from io import StringIO
import time

from Bio import SeqIO
import requests


def get_pLDDT(uid: str) -> list:
    '''
    Download pLDDT from AlphaFold.
    '''
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-confidence_v6.json'
    response = requests.get(url = url)

    if response.status_code != 200:
        return None
    return response.json()['confidenceScore']
    
def get_protein_sequences(uids: list, chunk_size: int = 500) -> str:
    '''
    Download protein FASTA sequence from UniPrtot.
    '''
    url = f"https://rest.uniprot.org/uniprotkb/accessions"
    result = {}
    for i in range(0, len(uids), 500):
        chunk = uids[i:i + chunk_size]
        response = requests.get(url = url, params = {'accessions': ','.join(chunk), 'format': 'tsv', 'fields': 'gene_primary, sequence'})
        for line in response.text.splitlines()[1:]:
            gene, seq = line.split("\t", 1)
            result[gene] = seq
    return result