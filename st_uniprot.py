import requests
import json
import streamlit as st

def get_human_uniprot_id(gene_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606&fields=accession"
    response = requests.get(url)
    data = response.json()
    return data['results'][0]['primaryAccession'] if data['results'] else None

def get_uniprot_length(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=length"
    response = requests.get(url)
    response.raise_for_status()
    data = json.loads(response.text)
    if data['results']:
        return data['results'][0]['sequence']['length']
    return "Length not found"

def get_uniprot_disease(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=cc_disease"
    response = requests.get(url)
    response.raise_for_status()
    data = json.loads(response.text)
    
    diseases = []
    if data['results'] and 'comments' in data['results'][0]:
        disease_comments = [comment for comment in data['results'][0]['comments'] if comment['commentType'] == 'DISEASE']
        for comment in disease_comments:
            disease = comment['disease']
            diseases.append(f"**_{disease['diseaseId']}_**. {disease['description']}")
    
    if diseases:
        return "\n".join(diseases)
    return "No disease information found"