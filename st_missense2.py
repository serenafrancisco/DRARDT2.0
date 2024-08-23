import subprocess
import tempfile
import os
import csv
import requests
from freesasa import *
import pandas as pd
import streamlit as st

def load_aa_properties(tsv_file):
    volume_dict = {}
    polarity_dict = {}

    # Read the TSV file
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            aa = row['AA']
            volume_dict[aa] = float(row['V'])
            polarity_dict[aa] = float(row['H'])
    
    return volume_dict, polarity_dict

def run_freesasa(input_pdb):
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_pdb_file:
        temp_pdb_file.write(input_pdb.read())
        temp_pdb_file.flush()  # Ensure all data is written to the file
        
        command = [
            "freesasa",
            "-L",
            "-n", "20",
            "-t", "8",
            "--radii", "naccess",
            temp_pdb_file.name
        ]
        
        try:
            # Execute the FreeSASA command and capture the output
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            output = result.stdout  # This captures the stdout as a string
            return output
        
        except subprocess.CalledProcessError as e:
            st.error(f"An error occurred while running FreeSASA: {e}")
        except FileNotFoundError:
            st.error("FreeSASA is not installed or not found in your PATH.")
    
    return None

def parse_freesasa_output(freesasa_output):
    data = []
    lines = freesasa_output.splitlines()
    
    # Extract relevant lines starting with "RES"
    for line in lines:
        if line.startswith("RES"):
            columns = line.split()
            if len(columns) >= 11:
                data.append({
                    "Res": columns[2],
                    "Chain": columns[1],
                    "ResNum": int(columns[3]),
                    "Area": float(columns[4]),
                    "RSA": float(columns[5]),
                    "Sidechain": float(columns[6]),
                    "Mainchain": float(columns[7]),
                    "NonPolar": float(columns[8]),
                    "AllPolar": float(columns[9]),
                    "RelTotal": float(columns[10])
                })
    
    # Convert the data to a pandas DataFrame
    df = pd.DataFrame(data)
    return df

def get_rsa_for_residue(freesasa_output, residue_number):
    """
    Extract the RSA value for a specific residue based on residue number from the FreeSASA output string.

    Parameters:
    - freesasa_output (str): The output string from the FreeSASA command.
    - residue_number (int): The residue number (e.g., 185).

    Returns:
    - rsa_value (float): The relative solvent accessible surface area of the residue, or None if not found.
    """
    rsa_value = None
    
    if freesasa_output:
        lines = freesasa_output.splitlines()
        for line in lines:
            # Skip comment lines and lines that don't match the expected format
            if line.startswith("REM") or not line.strip() or not line.startswith('RES'):
                continue
            
            # Split the line into columns
            columns = line.split()
            
            # Check if the line has enough columns and starts with 'RES'
            if len(columns) >= 11 and columns[0] == 'RES':
                try:
                    res_number = int(columns[3])  # Ensure no extra spaces when parsing the residue number
                except ValueError:
                    continue
                
                if res_number == residue_number:
                    try:
                        # Update the column index to where "All-atoms REL" is located in your specific output
                        rsa_value = float(columns[5])  # Ensure no extra spaces when parsing the RSA value
                    except ValueError:
                        rsa_value = None
                    break
    
    if rsa_value is None:
        print(f"\033[1mNo RSA value found for residue number {residue_number}.\033[0m")
    else:
        print(f"\033[1mRSA value for residue number {residue_number}:\033[0m {rsa_value}")
        if rsa_value > 20:
            print(f"Residue {residue_number} is exposed to the solvent.")
        else:
            print(f"Residue {residue_number} is buried.")
    
    return rsa_value

def calculate_ddG(mutation, RSA, volume_dict, polarity_dict):

    aa_wt = mutation[0]  # Wild-type amino acid
    aa_mut = mutation[-1]  # Mutated amino acid

    # Calculate Vdiff and Hdiff
    Vdiff = volume_dict[aa_mut] - volume_dict[aa_wt]
    Hdiff = polarity_dict[aa_mut] - polarity_dict[aa_wt]

    # Calculate SimBa-NI ΔΔG
    ddG = -1.64 + 1.9 * (RSA / 100) + 0.49 * (Vdiff/100) - 0.12 * Hdiff

    return ddG

def check_pdb_coverage(uniprot_id, pos_wt):
    # Build the URL for the UniProt API request
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=xref_pdb"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()  # Parse JSON response
    
    results = []  # Store PDB structures and their coverage
    position_covered = False  # Flag to check if the position is covered
    covered_structure = None  # To store the first structure that covers the position
    
    # Check if results are present and the PDB cross-references exist
    if data.get('results') and 'uniProtKBCrossReferences' in data['results'][0]:
        pdb_structures = data['results'][0]['uniProtKBCrossReferences']
        
        for pdb in pdb_structures:
            # Check if the cross-reference is for a PDB structure
            if pdb['database'] == 'PDB':
                structure = f"PDB ID: {pdb['id']}"
                
                # Extract coverage information
                coverage_info = pdb['properties'][2]['value']  # Example format: 'A=160-453,B=200-250'
                chains = coverage_info.split(',')
                
                # Check if the position falls within the chain coverage
                for chain in chains:
                    chain_parts = chain.split('=')
                    if len(chain_parts) == 2:
                        chain_range = chain_parts[1]
                        if '-' in chain_range:
                            # Handle range, e.g., '160-453'
                            start, end = map(int, chain_range.split('-'))
                            if start <= pos_wt <= end:
                                position_covered = True
                                covered_structure = structure  # Store the structure that covers the position
                                break
                
                results.append((structure, coverage_info))  # Append PDB structure and coverage info
                
            if position_covered:
                break  # Exit loop early if position is covered
    
    # Handle the case where no PDB cross-references exist
    else:
        return f"No PDB cross-references found for UniProt ID: {uniprot_id}", None
    
    return position_covered, covered_structure
