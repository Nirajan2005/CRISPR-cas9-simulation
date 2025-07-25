import random

# --- Core CRISPR functions (from your code, unchanged) ---

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

def count_mismatches(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))

def simulate_dna_repair(edited_dna, repair_mode="NHEJ", donor_template=None, cut_site=None, edit_length=3):
    if repair_mode == "NHEJ":
        # NHEJ: random small indel at cut site
        indel_type = random.choice(["del", "ins", "none"])
        if indel_type == "del":
            del_len = random.randint(1, min(3, len(edited_dna)-cut_site))
            return edited_dna[:cut_site] + edited_dna[cut_site+del_len:]
        elif indel_type == "ins":
            ins_seq = ''.join(random.choices("ATCG", k=random.randint(1,3)))
            return edited_dna[:cut_site] + ins_seq + edited_dna[cut_site:]
        else:
            return edited_dna
    elif repair_mode == "HDR" and donor_template:
        # HDR: replace region with donor template at cut site
        return edited_dna[:cut_site] + donor_template + edited_dna[cut_site+edit_length:]
    else:
        return edited_dna

def simulate_crispr_edit(
    dna, guide_rna, pam="NGG", edit="deletion", edit_length=3,
    insert_seq=None, substitute_seq=None, max_edits=None, strand="both",
    off_target_prob=0.1, max_mismatches=2, repair_mode="NHEJ", donor_template=None
):
    pam_sites = find_pam_sites(dna, pam)
    edits = []
    # Forward strand
    if strand in ("both", "forward"):
        for site in pam_sites:
            cut_site = site - 3
            target_start = cut_site - len(guide_rna) + 1
            if target_start < 0:
                continue
            target_seq = dna[target_start:cut_site+1]
            mismatches = count_mismatches(target_seq, guide_rna)
            is_on_target = mismatches == 0
            is_off_target = 0 < mismatches <= max_mismatches and random.random() < off_target_prob
            if is_on_target or is_off_target:
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
                # DNA repair outcome
                repaired_dna = simulate_dna_repair(
                    edited_dna, repair_mode=repair_mode, donor_template=donor_template,
                    cut_site=target_start + (0 if edit == "insertion" else 0), edit_length=edit_length
                )
                edits.append((repaired_dna, target_start, site, "forward", mismatches))
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
            target_seq = rev_dna[target_start:cut_site+1]
            mismatches = count_mismatches(target_seq, rev_guide)
            is_on_target = mismatches == 0
            is_off_target = 0 < mismatches <= max_mismatches and random.random() < off_target_prob
            if is_on_target or is_off_target:
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
                repaired_dna = simulate_dna_repair(
                    edited_dna, repair_mode=repair_mode, donor_template=donor_template,
                    cut_site=orig_target_start + (0 if edit == "insertion" else 0), edit_length=edit_length
                )
                edits.append((repaired_dna, orig_target_start, orig_cut_site, "reverse", mismatches))
                if max_edits and len(edits) >= max_edits:
                    return edits
    return edits

# --- Organ and DeliveryMethod classes ---

class DeliveryMethod:
    def __init__(self, name, organ_success_prob):
        self.name = name
        self.organ_success_prob = organ_success_prob  # dict: organ_name -> probability

    def get_success_prob(self, organ_name):
        return self.organ_success_prob.get(organ_name, 0.05)  # default low

class Organ:
    def __init__(self, name, dna_sequence, immune_clearance=0.1):
        self.name = name
        self.dna_sequence = dna_sequence
        self.immune_clearance = immune_clearance  # probability that delivery fails due to immune system

    def attempt_edit(self, delivery_method, **edit_kwargs):
        # Simulate immune clearance
        if random.random() < self.immune_clearance:
            return False, "Immune clearance", None
        # Simulate delivery success
        success_prob = delivery_method.get_success_prob(self.name)
        if random.random() > success_prob:
            return False, "Delivery failed", None
        # Attempt CRISPR edit
        edits = simulate_crispr_edit(self.dna_sequence, **edit_kwargs)
        if edits:
            return True, "Edit successful", edits[0]
        else:
            return False, "No target found", None

# --- Example setup ---

organs = [
    Organ("Liver", "ATGCGTACGTAGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGC", immune_clearance=0.05),
    Organ("Lung",  "CGTACGTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTA", immune_clearance=0.15),
    Organ("Muscle","GCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCT", immune_clearance=0.10),
    Organ("Brain", "TAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAGCTAGCTAGGCTAG", immune_clearance=0.20),
]

delivery_methods = [
    DeliveryMethod("LNP", {"Liver": 0.8, "Lung": 0.3, "Muscle": 0.2, "Brain": 0.05}),
    DeliveryMethod("AAV", {"Liver": 0.4, "Lung": 0.5, "Muscle": 0.5, "Brain": 0.7}),
    DeliveryMethod("Lentivirus", {"Liver": 0.3, "Lung": 0.4, "Muscle": 0.6, "Brain": 0.6}),
]

# --- User input ---

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

# --- Off-target and repair parameters ---
off_target_prob = float(input("Off-target probability (0-1, e.g. 0.1): ") or "0.1")
max_mismatches = int(input("Max mismatches allowed for off-target (e.g. 2): ") or "2")
repair_mode = input("DNA repair mode (NHEJ/HDR): ").upper() or "NHEJ"
donor_template = None
if repair_mode == "HDR":
    donor_template = input("Enter donor template for HDR (required): ").upper()

# --- Simulation ---

print("\n--- CRISPR Delivery Simulation Across Organs ---")
for organ in organs:
    print(f"\nOrgan: {organ.name}")
    for method in delivery_methods:
        print(f"  Delivery Method: {method.name}")
        # Safeguard for edit_length
        actual_edit_length = min(edit_length, len(organ.dna_sequence))
        success, reason, edit_result = organ.attempt_edit(
            method,
            guide_rna=guide_rna,
            edit=edit,
            edit_length=actual_edit_length,
            insert_seq=insert_seq,
            substitute_seq=substitute_seq,
            max_edits=max_edits,
            strand=strand,
            off_target_prob=off_target_prob,
            max_mismatches=max_mismatches,
            repair_mode=repair_mode,
            donor_template=donor_template
        )
        if success:
            edited_dna, cut_start, pam_site, strand_used, mismatches = edit_result
            if mismatches == 0:
                target_type = "ON-target"
            else:
                target_type = f"OFF-target ({mismatches} mismatch{'es' if mismatches>1 else ''})"
            print(f"    SUCCESS: {reason} [{target_type}]")
            print(f"    Edited DNA: {edited_dna}")
            print(f"    Edit at position {cut_start}-{pam_site} ({strand_used} strand)")
        else:
            print(f"    FAILURE: {reason}")