from Bio import Entrez
import requests
import json
import streamlit as st

def get_publication_count(gene_name):
    search_term = f"{gene_name} [Title/Abstract] AND 2000:2024 [PDat]"
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=1)
    record = Entrez.read(handle)
    return int(record["Count"])

def get_publication_count_score(pub_count):
    if pub_count <= 50:
        pub_count_score = 1
    elif pub_count <= 100:
        pub_count_score = 2
    elif pub_count <= 200:
        pub_count_score = 3
    else:
        pub_count_score = 4 
    return(pub_count_score)   

def get_string_interactors(gene_name):
    url = f"https://string-db.org/api/tsv-no-header/interaction_partners?identifiers={gene_name}&species=9606&network_type=physical"
    response = requests.get(url)
    response.raise_for_status()
    
    pathways = response.text.strip().split('\n')
    filtered_interactors = set()  # Use a set to avoid duplicates
    
    for pathway in pathways:
        if pathway:  # Ensure we don't process empty lines
            data = pathway.split('\t')
            experimental_score = float(data[10])  # 11th column (index 10) is the experimental score
            if experimental_score > 0.700:
                # Extract all protein names (columns 3 onward)
                for protein in data[2:4]:
                    if protein != gene_name:
                        filtered_interactors.add(protein)
    
    return len(filtered_interactors), list(filtered_interactors)

def get_interactors_score(interactors_count):
    if interactors_count < 2:
        interactors_score = 1
    elif interactors_count < 4:
        interactors_score = 2
    elif interactors_count <= 6:
        interactors_score = 3
    else:
        interactors_score = 4  
    return(interactors_score)  

def get_kegg_pathways(gene_name):
    url = f"https://rest.kegg.jp/find/pathway/{gene_name}"
    response = requests.get(url)
    response.raise_for_status()
    pathways = response.text.split('\n')
    return len(pathways) - 1  # subtracting 1 to account for empty string at end

def get_KEGG_score(pathways_count):
    if pathways_count < 1:
        KEGG_score = 1
    elif pathways_count < 2:
        KEGG_score = 2
    else:
        KEGG_score = 3 
    return(KEGG_score)   

def get_uniprot_3d(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=xref_pdb"
    response = requests.get(url)
    response.raise_for_status()
    data = json.loads(response.text)
    structures = []
    
    pdb_count = 0
    
    if data['results'] and 'uniProtKBCrossReferences' in data['results'][0]:
        pdb_structures = data['results'][0]['uniProtKBCrossReferences']
        for pdb in pdb_structures:
            if pdb['database'] == 'PDB':
                pdb_count += 1
                structure_info = (
                    f"PDB ID: {pdb['id']}, Method: {pdb['properties'][0]['value']}, "
                    f"Resolution: {pdb['properties'][1]['value']}, Chains: {pdb['properties'][2]['value']}"
                )
                structures.append(structure_info)
    
    if pdb_count == 0:
        structures.append("_No experimental 3D structures found_")

    PDB_score = 0
    if pdb_count < 1:
        PDB_score = 1
    elif pdb_count <= 3:
        PDB_score = 2
    elif pdb_count <= 4:
        PDB_score = 3
    else:
        PDB_score = 4
    
    return pdb_count, structures, PDB_score

def get_alphafold_prediction(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
    else:
        return "No AlphaFold prediction found"

def get_AF2_score(alphafold2_prediction):
    if alphafold2_prediction is not None:
        AF2_score = 2
    else:
        AF2_score = 1
    return(AF2_score)

def calculate_DRARDT_score(pub_count_score, interactors_score, KEGG_score, PDB_score, AF2_score):
    sum_scores = pub_count_score + interactors_score + KEGG_score + PDB_score + AF2_score
    if sum_scores < 8:
        DRARDT_score = 0
    elif 8 <= sum_scores < 11:
        DRARDT_score = 1
    elif 11 <= sum_scores < 14:
        DRARDT_score = 2
    elif sum_scores >= 14:
        DRARDT_score = 3
    return(DRARDT_score)