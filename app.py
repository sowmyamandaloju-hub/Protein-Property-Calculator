import streamlit as st
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

st.set_page_config(page_title="Protein Property Calculator", layout="centered")

st.title("🔬 Protein Property Calculator")
st.write("Paste a protein sequence (plain or FASTA) and get key physicochemical properties.")

def clean_sequence(seq):
    seq = seq.strip()

    # Remove FASTA header if present
    if seq.startswith(">"):
        lines = seq.split("\n")
        seq = "".join(lines[1:])

    seq = seq.replace("\n", "").replace(" ", "").upper()
    return seq

def calculate_properties(seq):
    analysis = ProteinAnalysis(seq)

    props = {
        "Sequence Length": len(seq),
        "Molecular Weight": round(analysis.molecular_weight(), 2),
        "Isoelectric Point (pI)": round(analysis.isoelectric_point(), 2),
        "Aromaticity": round(analysis.aromaticity(), 3),
        "Instability Index": round(analysis.instability_index(), 2),
        "GRAVY (Hydrophobicity)": round(analysis.gravy(), 3),
        "Extinction Coefficient (Reduced)": analysis.molar_extinction_coefficient()[0],
        "Extinction Coefficient (Oxidized)": analysis.molar_extinction_coefficient()[1],
    }

    aa_comp = analysis.get_amino_acids_percent()
    aa_df = pd.DataFrame(
        [{"Amino Acid": aa, "Percentage (%)": round(p * 100, 2)} for aa, p in aa_comp.items()]
    )

    props_df = pd.DataFrame([props])

    return props_df, aa_df


seq_input = st.text_area("🧬 Enter Protein Sequence", height=200)

if st.button("Calculate Properties"):
    seq = clean_sequence(seq_input)

    if seq == "":
        st.error("Please enter a protein sequence.")
    elif not seq.isalpha():
        st.error("Invalid sequence! Only amino acid letters allowed.")
    else:
        props_df, aa_df = calculate_properties(seq)

        st.success("✅ Calculation completed!")

        st.subheader("📌 Protein Properties")
        st.dataframe(props_df, use_container_width=True)

        st.subheader("🧬 Amino Acid Composition")
        st.dataframe(aa_df, use_container_width=True)

        # Download buttons
        st.download_button(
            "⬇️ Download Protein Properties CSV",
            props_df.to_csv(index=False),
            file_name="protein_properties.csv",
            mime="text/csv"
        )

        st.download_button(
            "⬇️ Download Amino Acid Composition CSV",
            aa_df.to_csv(index=False),
            file_name="amino_acid_composition.csv",
            mime="text/csv"
        )
