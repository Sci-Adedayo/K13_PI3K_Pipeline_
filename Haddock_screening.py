import os
import glob
import subprocess

def sanitize_pdb_chains():
    """Forces internal Chain B on targets and Chain A on mutants to prevent collisions."""
    print("[*] Sanitizing internal PDB chains...")
    
    # 1. Format Targets -> Chain B
    os.system("pdb_chain -B targets/PI3K_Catalytic_B.pdb > targets/tmp_cat.pdb && mv targets/tmp_cat.pdb targets/PI3K_Catalytic_B.pdb")
    os.system("pdb_chain -B targets/PI3K_Helical_B.pdb > targets/tmp_hel.pdb && mv targets/tmp_hel.pdb targets/PI3K_Helical_B.pdb")
    
    # 2. Format Mutants -> Chain A
    for mutant in glob.glob("mutants/*_A.pdb"):
        os.system(f'pdb_chain -A "{mutant}" > "{mutant}.tmp" && mv "{mutant}.tmp" "{mutant}"')
        
    print("[✓] PDB chain IDs successfully formatted.")

def run_local_haddock(mutant_pdb, domain, target_pdb, restraints_tbl, cores=16):
    mutant_name = os.path.basename(mutant_pdb).replace("_A.pdb", "")
    run_dir = f"run_{mutant_name}_{domain}"
    
    os.makedirs(run_dir, exist_ok=True)
    os.system(f"cp {mutant_pdb} {run_dir}/Receptor_A.pdb")
    os.system(f"cp {target_pdb} {run_dir}/Ligand_B.pdb")
    os.system(f"cp {restraints_tbl} {run_dir}/restraints.tbl")
    
    config = f"""run_dir = "output"
mode = "local"
ncores = {cores}

molecules = [
    "Receptor_A.pdb",
    "Ligand_B.pdb"
]

[topoaa]

[rigidbody]
sampling = 1000
ambig_fname = "restraints.tbl"

[seletop]
select = 200

[flexref]

[emref]

[caprieval]

[clustfcc]"""

    with open(f"{run_dir}/docking.cfg", "w") as f:
        f.write(config)
        
    print(f"[*] Starting {mutant_name} against {domain} domain...")
    
    cwd = os.getcwd()
    os.chdir(run_dir)
    subprocess.run(["haddock3", "docking.cfg"])
    os.chdir(cwd)

if __name__ == "__main__":
    mutant_files = sorted(glob.glob("mutants/*_A.pdb"))
    if not mutant_files:
        print("ERROR: No mutant structures found in mutants/ directory.")
        exit(1)
    
    # Automatically execute the PDB chain cleanup before docking
    sanitize_pdb_chains()
        
    print(f"Found {len(mutant_files)} mutants. Starting local batch screening...")
    
    for mutant in mutant_files:
        run_local_haddock(mutant, "catalytic", "targets/PI3K_Catalytic_B.pdb", "restraints/restraints_catalytic.tbl")
        run_local_haddock(mutant, "helical", "targets/PI3K_Helical_B.pdb", "restraints/restraints_helical.tbl")
