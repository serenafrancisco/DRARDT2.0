import streamlit as st
from Bio import Entrez
from st_uniprot import *
from st_params import *
from st_missense import *

def main():

    
    st.title("DRARDT.py")

    st.subheader("Web app of the _Drug Repurposing Assessment for Rare Disease Targets_ method")
    st.markdown("")

    cassmedchem_url = "https://www.cassmedchem.unito.it/"
    cassmedchem_logo = "cassmedchem_logo.png"
    paper_drardt_url = ""
    drardt_scheme = "drardt.jpg"

    st.markdown("DRARDT.py was developed at [**CASSMedChem**](%s) (University of Torino, Italy)." % cassmedchem_url)
    st.image(cassmedchem_logo)
    st.markdown("The original DRARDT method has been presented for the first time in [cite paper]:")
    st.image(drardt_scheme, caption="Adapted from [add reff paper]")
    st.link_button("Check out the original DRARDT publication", paper_drardt_url)

    st.markdown("\n\n")
    st.markdown("**DRARDT.py** requires the user to provide a gene name and email address (mandatory).")
    
    st.markdown(
    """
    With respect to the original version of the method, **DRARDT.py** takes into account the following aspects related to the target gene, and assigns an individual score to each one:
    - Number of publications in PubMed featuring the target from 2000 until present (and related score)
    - Number and names of interactors of the target protein from STRING-DB (and related score)
    - Number of KEGG pathways associated with the query gene (and related score)
    - Number and details of PDB entries for the target protein (and related score)
    - Availability of an AlphaFold2 prediction of the target protein (and related score)
    """
    )

    st.markdown("Additionally, this web app also reports information about the target UniProt ID, protein length and related disease(s) description.")

    st.markdown("Finally, the app prints a comprehensive (DRARDT) score from individual scores collected previously. "
                "This score can either be :red[**0** (Very Low)], " 
                ":orange[**1** (Low)], "
                ":green[**2** (High)], "
                "or :blue[**3** (Very High)]. "
                "Increasing scores indicate targets with higher potential of being targeted by small molecule-based drug repurposing.")
    
    drardt_expl = ""
    
    st.markdown("Details about the DRARDT scoring method are disclosed at this [link](%s)." % drardt_expl)


    gene_name = st.text_input("Enter gene name:")
    email = st.text_input("Enter a valid email address:")

    st.markdown("_OPTIONAL_. The user can also upload a list of missense mutations and a PDB file of the target "
                "to be checked for coverage in available PDB structures, solvent exposure (calculated with "
                "FreeSASA) and entity of downstream destabilization effect at the protein structure level.")
    mutations = st.text_input("Enter one or more missense mutations in the A123B format, separated by commas (e.g. A123B,C456D):")

    if mutations:
        input_pdb = st.file_uploader("Upload PDB file for the target (mandatory whenever missense mutations are submitted)", type="pdb")  # PDB required only when mutations are submitted

    if st.button("Submit"):
        if not email or not gene_name:
            st.write("Please provide both gene name and email address.")
        else:
            Entrez.email = email
            
            
            st.subheader("Generalities")
            
            uniprot_id = get_human_uniprot_id(gene_name)
            st.write(f":blue[**UNIPROT ID FOR {gene_name}:**]")
            st.write(uniprot_id)

            if uniprot_id:

                uniprot_pathology = get_uniprot_disease(uniprot_id)
                st.write(f":blue[**INVOLVEMENT OF {gene_name} IN DISEASES:**]")
                diseases = uniprot_pathology.split("\n")
                for disease in diseases:
                    st.write(disease)

                uniprot_prot_length = get_uniprot_length(uniprot_id)
                st.write(f":blue[**LENGTH OF {gene_name}:**]") 
                st.write(f"{uniprot_prot_length}")

                st.subheader("DRARDT parameters assessing")

                pub_count = get_publication_count(gene_name)
                st.write(f":blue[**NUMBER OF PUBLICATIONS ABOUT {gene_name} FROM 2000 T0 2024:**]") 
                st.write(f"{pub_count}")
                pub_count_score = get_publication_count_score(pub_count)
                st.write(f":violet[**Publication count score for {gene_name}:**]", (pub_count_score)) 

                interactors_count, interactors_names = get_string_interactors(gene_name)
                st.write(f":blue[**NUMBER OF {gene_name} INTERACTORS FROM STRING-DB:**]") 
                st.write(f"{interactors_count}")
                if interactors_count == 0:
                    st.write(f":blue[**INTERACTORS OF {gene_name}:**]") 
                    st.write("_None_")
                else:
                    st.write(f":blue[**INTERACTORS OF {gene_name}:**]") 
                    st.write(f"{', '.join(interactors_names)}")

                interactors_score = get_interactors_score(interactors_count)
                st.write(f":violet[**Interactors score for {gene_name}:**]", (interactors_score)) 

                pathways_count = get_kegg_pathways(gene_name)
                st.write(f":blue[**NUMBER OF KEGG PATHWAYS {gene_name} IS INVOLVED IN:**]") 
                st.write(f"{pathways_count}")
                KEGG_score = get_KEGG_score(pathways_count)
                st.write(f":violet[**KEGG score for {gene_name}:**]", (KEGG_score))

                pdb_count, structures, PDB_score = get_uniprot_3d(uniprot_id)
                st.write(f":blue[**NUMBER OF PDB ENTRIES FOR {gene_name}:**]") 
                st.write(f"{pdb_count}")
                for structure in structures:
                    st.write(structure)
                st.write(f":violet[**PDB score for {gene_name}:**]", (PDB_score))

                alphafold2_prediction = get_alphafold_prediction(uniprot_id)
                st.write(f":blue[**ALPHAFOLD2 PREDICTION FOR {gene_name} AVAILABLE AT:**]") 
                st.write(alphafold2_prediction)
                AF2_score = get_AF2_score(alphafold2_prediction)
                st.write(f":violet[**AlphaFold2 score for {gene_name}:** ]",(AF2_score))

                st.subheader(f"DRARDT score for {gene_name}")

                DRARDT_score = calculate_DRARDT_score(pub_count_score, interactors_score, KEGG_score, PDB_score, AF2_score)
                st.write(f"**Score:** {DRARDT_score}")
                if DRARDT_score == 0:
                    st.write(f"Flag: red:[**Very Low**]")
                elif DRARDT_score == 1:
                    st.write(f"Flag: :orange[**Low**]")   
                elif DRARDT_score == 2:
                    st.write(f"Flag: :green[**High**]")  
                elif DRARDT_score == 3:
                    st.write(f"Flag: :blue[**Very High**]\n")                

        if mutations and input_pdb:

            tsv_file = "simba.tsv"
            volume_dict, polarity_dict = load_aa_properties(tsv_file)

            mutations = [mutation.strip() for mutation in mutations.split(',')]
            for mutation in mutations:
                pos_wt = int(mutation[1:-1])
                is_covered, pdb_structure = check_pdb_coverage(uniprot_id, pos_wt)
                if is_covered:
                    st.write(f"Position {pos_wt} is covered by at least one experimental structure of {gene_name} ({pdb_structure}).")
                else:
                    st.write(f"Position {pos_wt} is not covered by any experimental structures of {gene_name}.")

            sasa_out = run_freesasa(gene_name, input_pdb)
            st.write(f"FreeSASA output saved to {sasa_out}")

            for mutation in mutations:

                pos_wt = int(mutation[1:-1])
                rsa_value = get_rsa_for_residue(sasa_out, pos_wt)

                if rsa_value is not None:
                    ddG = calculate_ddG(mutation, rsa_value, volume_dict, polarity_dict)
                    st.write(f"**Calculated SimBa-NI ΔΔG for {mutation}:** {ddG}")
                    if ddG > -1.5:
                        st.write(f"Mutation {mutation} is not expected to lead to protein unfolding.")
                    else:
                        st.write(f"Mutation {mutation} is expected to lead to protein unfolding.")
                else:
                    st.write(f"Could not calculate SimBa-NI ΔΔG for {mutation} due to missing RSA value.")

if __name__ == "__main__":
    main()