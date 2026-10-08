from io import StringIO

from Bio import SeqIO
import requests


def get_pLDDT(uid: str) -> list:
    '''
    Download pLDDT from AlphaFold.
    '''
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-confidence_v6.json'
    response = requests.get(url = url)

    return response.json()['confidenceScore']

def get_protein_sequence(uid: str) -> str:
    '''
    Download protein FASTA sequence from UniPrtot.
    '''
    url = f'https://rest.uniprot.org/uniprotkb/{uid}.fasta'
    response = requests.get(url = url)
    
    fasta_io = StringIO(response.text)
    return str(SeqIO.read(fasta_io, 'fasta').seq)