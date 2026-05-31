"""Codon-aware mutation detection for DNA sequences relative to a wild-type reference."""

from typing import Dict, List, Tuple

# Standard genetic code mapping codons to single-letter amino acids
GENETIC_CODE: Dict[str, str] = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "_", "TAG": "_",
    "TGC": "C", "TGT": "C", "TGA": "_", "TGG": "W",
}

# Mapping of single-letter amino acids to standard three-letter abbreviations
AA_THREE_LETTER: Dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "_": "Stop", "?": "Unk"
}


def translate_codon(codon: str) -> str:
    """Translate a 3-nucleotide DNA codon to its 3-letter amino acid code.

    If the codon has gaps or is invalid, returns 'Frameshift' or 'Unk'.
    """
    cleaned = codon.replace("-", "").upper()
    if len(cleaned) != 3:
        return "Frameshift"
    aa_single = GENETIC_CODE.get(cleaned, "?")
    return AA_THREE_LETTER.get(aa_single, "Unk")


def detect_mutations(aligned_ref: str, aligned_mut: str) -> List[Dict]:
    """Scan aligned reference and mutated sequences to flag point mutations at the codon level.

    Args:
        aligned_ref: The aligned reference sequence (with gaps).
        aligned_mut: The aligned mutated sequence (with gaps).

    Returns:
        A list of dictionaries representing individual mutations.
    """
    # 1. Map each column of the alignment to a reference codon index (1-based)
    alignment_len = len(aligned_ref)
    codon_map = [0] * alignment_len
    ref_char_count = 0

    # First pass: map columns containing reference nucleotides
    for col in range(alignment_len):
        if aligned_ref[col] != "-":
            codon_map[col] = (ref_char_count // 3) + 1
            ref_char_count += 1

    # Second pass: fill in gaps (insertions in mutant relative to reference)
    # Assign them to the adjacent reference codon
    last_codon = 1
    for col in range(alignment_len):
        if codon_map[col] > 0:
            last_codon = codon_map[col]
        else:
            codon_map[col] = last_codon

    # Determine maximum reference codon number
    max_codon = (ref_char_count + 2) // 3

    # Group alignment characters by reference codon
    codon_columns: Dict[int, List[Tuple[str, str]]] = {
        c: [] for c in range(1, max_codon + 1)
    }
    for col in range(alignment_len):
        c_idx = codon_map[col]
        if c_idx in codon_columns:
            codon_columns[c_idx].append((aligned_ref[col], aligned_mut[col]))

    mutations = []

    # 2. Analyze each codon group
    for codon_idx in range(1, max_codon + 1):
        cols = codon_columns[codon_idx]
        ref_seg = "".join(r for r, m in cols)
        mut_seg = "".join(m for r, m in cols)

        # Remove trailing/leading gaps in the reference codon segment to get original WT codon
        ref_codon = ref_seg.replace("-", "")
        mut_codon_clean = mut_seg.replace("-", "")

        # Skip if there's no change
        if ref_seg == mut_seg:
            continue

        ref_aa = translate_codon(ref_codon)

        # Case A: Substitutions (lengths match and no gaps)
        if "-" not in ref_seg and "-" not in mut_seg:
            mut_aa = translate_codon(mut_seg)
            if ref_aa == mut_aa:
                mut_type = "Synonymous Substitution"
            elif mut_aa == "Stop":
                mut_type = "Nonsense Substitution"
            else:
                mut_type = "Missense Substitution"

            # Describe nucleotide substitution
            sub_details = []
            for i in range(len(ref_seg)):
                if ref_seg[i] != mut_seg[i]:
                    pos_in_codon = i + 1
                    sub_details.append(f"{ref_seg[i]}->{mut_seg[i]} at pos {pos_in_codon}")

            desc = f"p.{ref_aa}{codon_idx}{mut_aa} ({ref_seg}->{mut_seg})"
            mutations.append({
                "codon_number": codon_idx,
                "type": mut_type,
                "ref_seq": ref_seg,
                "mut_seq": mut_seg,
                "ref_amino_acid": ref_aa,
                "mut_amino_acid": mut_aa,
                "description": desc,
                "details": ", ".join(sub_details)
            })

        # Case B: Deletions (mutant segment contains gaps)
        elif "-" in mut_seg and "-" not in ref_seg:
            # Count deleted nucleotides
            deleted_nt = ref_seg.replace(mut_seg.replace("-", ""), "")
            # If the length of deleted characters is 0 (due to complex alignment), compute manually
            del_len = mut_seg.count("-")
            is_inframe = (del_len % 3 == 0)
            mut_type = "In-frame Deletion" if is_inframe else "Frameshift Deletion"
            mut_aa = "Deletion" if is_inframe else "Frameshift"

            desc = f"c.{codon_idx}_del{del_len} (Ref: {ref_seg} | Mut: {mut_seg})"
            mutations.append({
                "codon_number": codon_idx,
                "type": mut_type,
                "ref_seq": ref_seg,
                "mut_seq": mut_seg,
                "ref_amino_acid": ref_aa,
                "mut_amino_acid": mut_aa,
                "description": desc,
                "details": f"Deleted {del_len} bp ({ref_seg} -> {mut_seg})"
            })

        # Case C: Insertions (reference segment contains gaps)
        elif "-" in ref_seg and "-" not in mut_seg:
            inserted_nt = mut_seg.replace(ref_seg.replace("-", ""), "")
            ins_len = ref_seg.count("-")
            is_inframe = (ins_len % 3 == 0)
            mut_type = "In-frame Insertion" if is_inframe else "Frameshift Insertion"
            mut_aa = "Insertion" if is_inframe else "Frameshift"

            desc = f"c.{codon_idx}_ins{inserted_nt} (Ref: {ref_seg} | Mut: {mut_seg})"
            mutations.append({
                "codon_number": codon_idx,
                "type": mut_type,
                "ref_seq": ref_seg,
                "mut_seq": mut_seg,
                "ref_amino_acid": ref_aa,
                "mut_amino_acid": mut_aa,
                "description": desc,
                "details": f"Inserted '{inserted_nt}' ({ref_seg} -> {mut_seg})"
            })

        # Case D: Mixed / Complex Indels (both contain gaps)
        else:
            mut_type = "Complex Mutation (Indel)"
            mut_aa = "Complex"
            desc = f"c.{codon_idx}_complex (Ref: {ref_seg} -> Mut: {mut_seg})"
            mutations.append({
                "codon_number": codon_idx,
                "type": mut_type,
                "ref_seq": ref_seg,
                "mut_seq": mut_seg,
                "ref_amino_acid": ref_aa,
                "mut_amino_acid": mut_aa,
                "description": desc,
                "details": f"Ref aligned as '{ref_seg}', Mutant as '{mut_seg}'"
            })

    return mutations
