import random

def find_pam_sites(dna, pam="NGG"):
    pam_sites = []
    for i in range(len(dna) - len(pam) + 1):
        match = True
        for j, base in enumerate(pam):
            if base != "N" and dna[i + j] != base:
                match = False
                break
        if match:
            pam_sites.append(i)
    return pam_sites

def reverse_complement(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(comp[base] for base in reversed(seq))

def simulate_crispr_edit(
    dna, guide_rna, pam="NGG", edit="deletion", edit_length=3,
    insert_seq=None, substitute_seq=None, max_edits=None, strand="both"
):
    pam_sites = find_pam_sites(dna, pam)
    edits = []
    strands_checked = []
    # Forward strand
    if strand in ("both", "forward"):
        for site in pam_sites:
            cut_site = site - 3  # Cas9 cuts ~3 bases upstream of PAM
            target_start = cut_site - len(guide_rna) + 1
            if target_start < 0:
                continue
            if dna[target_start:cut_site+1] == guide_rna:
                if edit == "deletion":
                    edited_dna = dna[:target_start] + dna[target_start + edit_length:]
                elif edit == "insertion":
                    if not insert_seq:
                        insert_seq = ''.join(random.choices("ATCG", k=edit_length))
                    edited_dna = dna[:target_start] + insert_seq + dna[target_start:]
                elif edit == "substitution":
                    if not substitute_seq:
                        substitute_seq = ''.join(random.choices("ATCG", k=edit_length))
                    edited_dna = dna[:target_start] + substitute_seq + dna[target_start + edit_length:]
                else:
                    edited_dna = dna
                edits.append((edited_dna, target_start, site, "forward"))
                if max_edits and len(edits) >= max_edits:
                    return edits
    # Reverse strand
    if strand in ("both", "reverse"):
        rev_dna = reverse_complement(dna)
        rev_pam_sites = find_pam_sites(rev_dna, pam)
        rev_guide = reverse_complement(guide_rna)
        for site in rev_pam_sites:
            cut_site = site - 3
            target_start = cut_site - len(rev_guide) + 1
            if target_start < 0:
                continue
            if rev_dna[target_start:cut_site+1] == rev_guide:
                # Map back to original coordinates
                orig_target_start = len(dna) - (cut_site + 1)
                orig_cut_site = len(dna) - (site)
                if edit == "deletion":
                    edited_dna = dna[:orig_target_start] + dna[orig_target_start + edit_length:]
                elif edit == "insertion":
                    if not insert_seq:
                        insert_seq = ''.join(random.choices("ATCG", k=edit_length))
                    edited_dna = dna[:orig_target_start] + insert_seq + dna[orig_target_start:]
                elif edit == "substitution":
                    if not substitute_seq:
                        substitute_seq = ''.join(random.choices("ATCG", k=edit_length))
                    edited_dna = dna[:orig_target_start] + substitute_seq + dna[orig_target_start + edit_length:]
                else:
                    edited_dna = dna
                edits.append((edited_dna, orig_target_start, orig_cut_site, "reverse"))
                if max_edits and len(edits) >= max_edits:
                    return edits
    return edits

# User input section
original_dna = input("Enter the DNA sequence: ").upper()
guide_rna = input("Enter the guide RNA sequence: ").upper()
edit = input("Choose edit type (deletion/insertion/substitution): ").lower()
edit_length = int(input("Enter edit length: "))
strand = input("Strand to edit? (forward/reverse/both): ").lower()
max_edits = input("Max edits to perform (blank for all): ")
max_edits = int(max_edits) if max_edits.strip() else None

insert_seq = None
substitute_seq = None
if edit == "insertion":
    insert_seq = input("Enter sequence to insert (leave blank for random): ").upper() or None
elif edit == "substitution":
    substitute_seq = input("Enter sequence to substitute (leave blank for random): ").upper() or None

# Safeguard for edit_length
if edit_length > len(original_dna):
    print("Edit length is longer than DNA sequence. Adjusting to sequence length.")
    edit_length = len(original_dna)

edits = simulate_crispr_edit(
    original_dna, guide_rna, edit=edit, edit_length=edit_length,
    insert_seq=insert_seq, substitute_seq=substitute_seq,
    max_edits=max_edits, strand=strand
)

if not edits:
    print("No edit occurred.")
else:
    for idx, (edited_dna, cut_start, pam_site, strand_used) in enumerate(edits, 1):
        print(f"\nEdit #{idx} ({strand_used} strand):")
        print("Original DNA: ", original_dna)
        print("Edited DNA:   ", edited_dna)
        if cut_start is not None:
            region_start = max(0, cut_start-5)
            region_end = min(len(original_dna), pam_site+5)
            print("Edit target region:")
            print(original_dna[region_start:region_end])
            print(" "*(cut_start-region_start) + "^"*(pam_site-cut_start+1))
            print(f"Edit occurred at position {cut_start}-{pam_site}")
        else:
            print("No edit occurred.")