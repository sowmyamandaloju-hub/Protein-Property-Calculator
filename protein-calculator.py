from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

def clean_sequence(seq):
    seq = seq.strip()

    # Remove FASTA header if present
    if seq.startswith(">"):
        lines = seq.split("\n")
        seq = "".join(lines[1:])

    seq = seq.replace("\n", "").replace(" ", "").upper()
    return seq

def protein_properties(seq):
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
    return props, aa_comp

if __name__ == "__main__":
    print("\n🔬 Protein Property Calculator\n")

    seq_input = input("Enter protein sequence (plain or FASTA):\n")
    seq = clean_sequence(seq_input)

    if not seq.isalpha():
        print("\n❌ Invalid sequence! Only letters allowed.")
        exit()

    props, comp = protein_properties(seq)

    print("\n✅ Protein Properties:\n")
    for k, v in props.items():
        print(f"{k}: {v}")

    print("\n🧬 Amino Acid Composition (%):\n")
    for aa, percent in comp.items():
        print(f"{aa}: {round(percent*100, 2)}%")

    # Save outputs
    pd.DataFrame([props]).to_csv("protein_properties.csv", index=False)
    pd.DataFrame(
        [{"Amino Acid": aa, "Percentage": round(p*100, 2)} for aa, p in comp.items()]
    ).to_csv("amino_acid_composition.csv", index=False)

    print("\n📁 Saved outputs:")
    print(" - protein_properties.csv")
    print(" - amino_acid_composition.csv")
