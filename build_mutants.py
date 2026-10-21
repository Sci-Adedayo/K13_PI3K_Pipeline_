#!/usr/bin/env python3
"""
PfK13 Mutant Structure Generation Script
========================================
Downloads 4ZGC (PfK13 propeller domain), isolates Chain A, and uses PyMOL 
to build wild-type and mutant PDB files formatted for downstream HADDOCK3 docking.
Enforces validated control mutations (C580Y, G449A) alongside CSV variants.
"""

import os
import re
import urllib.request
import pandas as pd
from pymol import cmd

OUTPUT_DIR = "mutants"
CSV_FILE = "K13 mutations and prevalence.csv"
PDB_WT = "4zgc.pdb"

# Enforce validated control mutations (Artemisinin resistance markers)
CONTROL_MUTATIONS = ["C580Y", "G449A"]

os.makedirs(OUTPUT_DIR, exist_ok=True)

AA_MAP = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL',
    'W': 'TRP', 'Y': 'TYR'
}

# -----------------------------------------------------------------------------
# 0. Download 4ZGC PDB if not locally present
# -----------------------------------------------------------------------------
if not os.path.isfile(PDB_WT):
    print(f"[+] Downloading Wild-Type structure ({PDB_WT})...")
    urllib.request.urlretrieve("https://files.rcsb.org/download/4ZGC.pdb", PDB_WT)

# -----------------------------------------------------------------------------
# 1. Prepare and Save Cleaned Isolated Wild-Type Control (Chain A)
# -----------------------------------------------------------------------------
print("[+] Processing and saving Wild-Type Control (Chain A)...")
cmd.reinitialize()
cmd.load(PDB_WT, "wt_complex")
cmd.remove("not chain A")
cmd.remove("not polymer.protein")
cmd.remove("not alt '' and not alt 'A'")
cmd.alter("all", "alt=''")
cmd.alter("all", "chain='A'")

wt_out_path = os.path.join(OUTPUT_DIR, "K13_WT_A.pdb")
cmd.save(wt_out_path, "wt_complex")
print(f"[✓] Wild-Type control saved to: {wt_out_path}")

# -----------------------------------------------------------------------------
# 2. Load Mutations from CSV & Append Validated Controls (C580Y, G449A)
# -----------------------------------------------------------------------------
raw_mutations = []

if os.path.isfile(CSV_FILE):
    try:
        df = pd.read_csv(CSV_FILE)
        if 'Mutation' in df.columns:
            raw_mutations = df['Mutation'].dropna().str.strip().unique().tolist()
            print(f"[+] Loaded {len(raw_mutations)} mutations from '{CSV_FILE}'.")
        else:
            print(f"[!] Warning: Column 'Mutation' not found in '{CSV_FILE}'.")
    except Exception as e:
        print(f"[!] Error reading '{CSV_FILE}': {e}")
else:
    print(f"[!] '{CSV_FILE}' not found locally. Generating control mutations only.")

# Merge CSV mutations with validated controls while preventing duplicates
all_mutations = list(raw_mutations)
for ctrl in CONTROL_MUTATIONS:
    if ctrl not in all_mutations:
        all_mutations.append(ctrl)
        print(f"[+] Enforced validated control mutation: {ctrl}")

all_mutations = sorted(list(set(all_mutations)))
print(f"[+] Total variant structures to generate (including controls): {len(all_mutations)}")

# -----------------------------------------------------------------------------
# 3. Build Mutant PDB Structures from Isolated Chain A
# -----------------------------------------------------------------------------
for mut in all_mutations:
    match = re.match(r"([A-Z])(\d+)([A-Z])", mut)
    if not match:
        print(f"[!] Skipping invalid mutation format: '{mut}'")
        continue
        
    wt_aa, res_num, target_aa = match.groups()
    res_num = int(res_num)
    target_aa_3 = AA_MAP.get(target_aa)
    
    if not target_aa_3:
        print(f"[!] Skipping unknown target amino acid in '{mut}'")
        continue
    
    out_path = os.path.join(OUTPUT_DIR, f"K13_{mut}_A.pdb")
    
    if os.path.isfile(out_path):
        print(f"[*] Skipping {mut}: Output structure already exists at {out_path}")
        continue

    print(f"[*] Building mutant {mut} (Residue {res_num} -> {target_aa_3})...")
    
    cmd.reinitialize()
    cmd.load(PDB_WT, "wt_complex")
    cmd.remove("not chain A")
    cmd.remove("not polymer.protein")
    cmd.remove("not alt '' and not alt 'A'")
    cmd.alter("all", "alt=''")
    cmd.alter("all", "chain='A'")
    
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    cmd.get_wizard().set_mode(target_aa_3)
    cmd.get_wizard().do_select(f"resi {res_num}")
    cmd.get_wizard().apply()
    cmd.set_wizard()
    
    cmd.save(out_path, "wt_complex")

cmd.quit()
print(f"\n[✓] Finished! Wild-Type control (K13_WT_A.pdb) and {len(all_mutations)} mutant structures are ready in '{OUTPUT_DIR}/'.")