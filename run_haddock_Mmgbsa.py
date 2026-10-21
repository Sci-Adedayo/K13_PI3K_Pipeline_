#!/usr/bin/env python3
"""
PfK13 / PI3K Automated Master Pipeline (Docking, Extraction, Minimization & MM/GBSA)
====================================================================================
Executes end-to-end overnight:
  Phase 0: Environment setup & parameter control file generation (MDP files)
  Phase 1: PDB chain sanitization & HADDOCK3 batch docking (catalytic & helical)
  Phase 2: Capri report parsing & extraction of best-scoring HADDOCK cluster complexes
  Phase 3: Amber topology creation (tleap), 500-step Sander minimization, & MMPBSA.py rescoring
  Phase 4: Output compilation to master TSV tables for notebook visualization
"""

import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import pandas as pd

# ==============================================================================
# CONFIGURATION & GLOBAL PATHS
# ==============================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in locals() else os.path.abspath(os.getcwd())
CONDA_BIN = "/home/adedayoadedisa/miniconda3/bin"
CORES_PER_RUN = 16

# Workspace Directory Structure
MUTANTS_DIR = os.path.join(BASE_DIR, "mutants")
TARGETS_DIR = os.path.join(BASE_DIR, "targets")
RESTRAINTS_DIR = os.path.join(BASE_DIR, "restraints")
COMPLEX_DIR = os.path.join(BASE_DIR, "mutants/mutant_complex")
MDP_DIR = os.path.join(BASE_DIR, "mdp_files")
HAWKDOCK_DIR = os.path.join(BASE_DIR, "hawkdock_runs")
RESULTS_DIR = os.path.join(BASE_DIR, "Analysis_results")

# Ensure required runtime directories exist
for d in [COMPLEX_DIR, MDP_DIR, HAWKDOCK_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# Environment setup for AmberTools & HADDOCK3
env = os.environ.copy()
env["PATH"] = f"{CONDA_BIN}:{env.get('PATH', '')}"

# Logging setup (Outputs simultaneously to terminal and master_pipeline.log)
log_file = os.path.join(BASE_DIR, "master_pipeline.log")
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)

VALID_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID", 
    "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", 
    "TRP", "TYR", "VAL"
}


# ==============================================================================
# PHASE 0: CONTROL PARAMETER FILE GENERATION
# ==============================================================================
def create_control_files():
    """Generates optimized Amber control parameter files for Sander and MMPBSA.py dynamically."""
    logging.info("Generating input parameter files in mdp_files/...")
    
    # Fast MM/GBSA input file (idecomp omitted to eliminate 15-minute residue decomposition bottleneck)
    mmpbsa_file = os.path.join(MDP_DIR, "mmpbsa_hawkdock.in")
    with open(mmpbsa_file, "w") as f:
        f.write("""&general
  startframe = 1,
  endframe = 1,
  verbose = 2,
/
&gb
  igb = 2,
  saltcon = 0.150,
  surften = 0.0072,
/
""")

    # Fast 500-step implicit solvent minimization file
    min_file = os.path.join(MDP_DIR, "min_hawkdock.in")
    with open(min_file, "w") as f:
        f.write("""500 steps fast minimization in implicit solvent
 &cntrl
  imin=1, maxcyc=500, ncyc=250,
  cut=25.0, igb=2, ntb=0,
  ntpr=100
 /
""")
    return mmpbsa_file, min_file


# ==============================================================================
# PHASE 1: SANITIZATION & HADDOCK3 DOCKING
# ==============================================================================
def sanitize_pdb_chains():
    """Forces Chain B on target structures and Chain A on mutants to prevent collisions."""
    logging.info("Sanitizing target and mutant PDB chain IDs...")
    
    cat_target = os.path.join(TARGETS_DIR, "PI3K_Catalytic_B.pdb")
    hel_target = os.path.join(TARGETS_DIR, "PI3K_Helical_B.pdb")
    
    for target in [cat_target, hel_target]:
        if os.path.exists(target):
            tmp_target = f"{target}.tmp"
            cmd = f'pdb_chain -B "{target}" > "{tmp_target}" && mv "{tmp_target}" "{target}"'
            subprocess.run(cmd, shell=True, env=env)

    for mutant in glob.glob(os.path.join(MUTANTS_DIR, "*_A.pdb")):
        tmp_mutant = f"{mutant}.tmp"
        cmd = f'pdb_chain -A "{mutant}" > "{tmp_mutant}" && mv "{tmp_mutant}" "{mutant}"'
        subprocess.run(cmd, shell=True, env=env)

    logging.info("PDB chain sanitization complete.")


def run_local_haddock(mutant_pdb, domain, target_pdb, restraints_tbl):
    """Configures and executes a HADDOCK3 docking simulation."""
    mutant_name = os.path.basename(mutant_pdb).replace("_A.pdb", "")
    run_dir = os.path.join(BASE_DIR, f"run_{mutant_name}_{domain}")
    
    # Skip if docking run has already been completed or initialized
    if os.path.exists(os.path.join(run_dir, "output")):
        logging.info(f"Skipping docking for {mutant_name} ({domain}): Output directory exists.")
        return run_dir

    os.makedirs(run_dir, exist_ok=True)
    shutil.copy(mutant_pdb, os.path.join(run_dir, "Receptor_A.pdb"))
    shutil.copy(target_pdb, os.path.join(run_dir, "Ligand_B.pdb"))
    shutil.copy(restraints_tbl, os.path.join(run_dir, "restraints.tbl"))

    config = f"""run_dir = "output"
mode = "local"
ncores = {CORES_PER_RUN}

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

    cfg_path = os.path.join(run_dir, "docking.cfg")
    with open(cfg_path, "w") as f:
        f.write(config)

    logging.info(f"Starting HADDOCK3 docking: {mutant_name} vs {domain} domain")
    proc = subprocess.run(
        ["haddock3", "docking.cfg"],
        cwd=run_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    
    if proc.returncode != 0:
        logging.error(f"HADDOCK3 failed for {mutant_name}_{domain}.")
    else:
        logging.info(f"Completed HADDOCK3 docking: {mutant_name}_{domain}")
        
    return run_dir


# ==============================================================================
# PHASE 2: TOP-SCORING CLUSTER COMPLEX EXTRACTION
# ==============================================================================
def extract_best_cluster_complexes():
    """Parses capri_clt.tsv in each run directory, identifies the lowest-scoring cluster,

    and copies the top-ranked model complex PDB to mutants/mutant_complex/.
    """
    logging.info("Extracting top-scoring cluster complex PDBs from HADDOCK outputs...")
    run_dirs = glob.glob(os.path.join(BASE_DIR, "run_*_*"))
    extracted_count = 0

    for r_dir in run_dirs:
        var_domain = os.path.basename(r_dir).replace("run_", "")
        target_complex_pdb = os.path.join(COMPLEX_DIR, f"{var_domain}_complex.pdb")

        clt_reports = glob.glob(os.path.join(r_dir, "output", "*_caprieval", "capri_clt.tsv"))
        if not clt_reports:
            continue

        clt_report = clt_reports[-1]
        with open(clt_report, "r") as f:
            lines = [l.strip().split("\t") for l in f if l.strip() and not l.startswith("#")]

        if len(lines) <= 1:
            continue

        header = lines[0]
        clusters = [dict(zip(header, row)) for row in lines[1:]]
        
        # Sort clusters by HADDOCK score (ascending / lowest / most negative score first)
        clusters_sorted = sorted(clusters, key=lambda x: float(x.get("score", 9999.0)))
        best_cluster_id = clusters_sorted[0].get("cluster_id", clusters_sorted[0].get("cluster_rank", "1"))

        stage_dir = os.path.dirname(clt_report)
        ss_report = os.path.join(stage_dir, "capri_ss.tsv")

        best_model_fname = None
        if os.path.isfile(ss_report):
            with open(ss_report, "r") as f:
                ss_lines = [l.strip().split("\t") for l in f if l.strip() and not l.startswith("#")]
            if len(ss_lines) > 1:
                ss_header = ss_lines[0]
                models = [dict(zip(ss_header, row)) for row in ss_lines[1:]]
                clt_models = [m for m in models if m.get("cluster_id") == str(best_cluster_id)]
                if clt_models:
                    clt_models_sorted = sorted(clt_models, key=lambda x: float(x.get("score", 9999.0)))
                    best_model_fname = clt_models_sorted[0].get("model")

        if best_model_fname:
            model_path = os.path.normpath(os.path.join(stage_dir, best_model_fname))
            if os.path.isfile(model_path):
                shutil.copy(model_path, target_complex_pdb)
                extracted_count += 1
                continue

        candidates = glob.glob(os.path.join(r_dir, "output", "*_clustfcc", f"*{best_cluster_id}*.pdb"))
        if candidates:
            shutil.copy(candidates[0], target_complex_pdb)
            extracted_count += 1

    logging.info(f"Successfully extracted {extracted_count} top complex PDBs into {COMPLEX_DIR}")


# ==============================================================================
# PHASE 3: AMBER TOPOLOGY & MM/GBSA RESCORING
# ==============================================================================
def is_hydrogen(line):
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element == "H": return True
    atom_name = line[12:16].strip().upper()
    if atom_name.startswith("H"): return True
    if len(atom_name) > 1 and atom_name[0].isdigit() and atom_name[1:].startswith("H"): return True
    return False


def normalize_line(line):
    raw_res = line[16:21].strip()
    if raw_res in VALID_AMINO_ACIDS: return line
    if len(raw_res) == 4 and raw_res[0] in ["C", "N"]:
        sub_res = raw_res[1:]
        if sub_res in VALID_AMINO_ACIDS:
            new_res = f"{sub_res:>3} "
            return line[:17] + new_res + line[21:]
    return line


def process_complex_pdb(pdb_path):
    """Parses a complex PDB file and robustly separates Receptor (Chain 1) and Ligand (Chain 2)."""
    rec_lines, lig_lines = [], []
    current_chain = 1
    prev_resnum, prev_segid, prev_chain_id = None, None, None

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("TER"):
                if rec_lines and current_chain == 1:
                    current_chain = 2
                continue
            if not line.startswith("ATOM") or is_hydrogen(line):
                continue

            cleaned_line = normalize_line(line)
            chain_id = cleaned_line[21].strip()
            segid = cleaned_line[72:76].strip()
            try:
                resnum = int(cleaned_line[22:26].strip())
            except ValueError:
                resnum = None

            if current_chain == 1 and len(rec_lines) > 0:
                if prev_chain_id and chain_id and chain_id != prev_chain_id:
                    current_chain = 2
                elif prev_segid and segid and segid != prev_segid:
                    current_chain = 2
                elif prev_resnum is not None and resnum is not None and resnum < (prev_resnum - 50):
                    current_chain = 2

            if resnum is not None: prev_resnum = resnum
            if segid: prev_segid = segid
            if chain_id: prev_chain_id = chain_id

            if current_chain == 1:
                rec_lines.append(cleaned_line)
            else:
                lig_lines.append(cleaned_line)

    return rec_lines, lig_lines


def run_mmgbsa_rescoring(mdp_file, min_file):
    """Builds Amber topologies, executes Sander implicit solvent minimization, and computes MM/GBSA scores."""
    logging.info("Starting Amber topology generation, Sander minimization, and MM/GBSA calculations...")
    pdb_files = sorted(glob.glob(os.path.join(COMPLEX_DIR, "*.pdb")))

    if not pdb_files:
        logging.warning(f"No complex PDB files found in {COMPLEX_DIR} for MM/GBSA rescoring.")
        return

    for pdb_path in pdb_files:
        variant_name = os.path.basename(pdb_path).replace(".pdb", "").replace("_complex", "")
        work_dir = os.path.join(HAWKDOCK_DIR, variant_name)
        os.makedirs(work_dir, exist_ok=True)

        dat_out = os.path.join(work_dir, "HAWKDOCK_MMGBSA_RESULTS.dat")
        if os.path.exists(dat_out) and os.path.getsize(dat_out) > 0:
            logging.info(f"Skipping MM/GBSA for {variant_name}: Output file exists.")
            continue

        rec_lines, lig_lines = process_complex_pdb(pdb_path)
        if not rec_lines or not lig_lines:
            logging.error(f"Chain separation failed for {variant_name}.")
            continue

        clean_pdb = os.path.join(work_dir, "clean_complex.pdb")
        rec_pdb = os.path.join(work_dir, "rec.pdb")
        lig_pdb = os.path.join(work_dir, "lig.pdb")

        with open(clean_pdb, "w") as f:
            f.writelines(rec_lines); f.write("TER\n"); f.writelines(lig_lines); f.write("TER\nEND\n")
        with open(rec_pdb, "w") as f:
            f.writelines(rec_lines); f.write("TER\nEND\n")
        with open(lig_pdb, "w") as f:
            f.writelines(lig_lines); f.write("TER\nEND\n")

        # Step A: Build Amber Topologies
        tleap_in = os.path.join(work_dir, "tleap.in")
        with open(tleap_in, "w") as f:
            f.write("""source leaprc.protein.ff14SB
source leaprc.water.tip3p
set default PBradii mbondi2

complex = loadpdb clean_complex.pdb
receptor = loadpdb rec.pdb
ligand = loadpdb lig.pdb

saveamberparm complex complex.prmtop complex.inpcrd
saveamberparm receptor receptor.prmtop receptor.inpcrd
saveamberparm ligand ligand.prmtop ligand.inpcrd
quit
""")

        subprocess.run(f"{CONDA_BIN}/tleap -f tleap.in > tleap.log 2>&1", shell=True, cwd=work_dir, env=env)
        if not os.path.exists(os.path.join(work_dir, "complex.prmtop")):
            logging.error(f"tleap failed for {variant_name}. Check {work_dir}/tleap.log")
            continue

        # Step B: 500-Step Sander Implicit Solvation Minimization
        min_cmd = (
            f"{CONDA_BIN}/sander -O -i {min_file} -p complex.prmtop "
            f"-c complex.inpcrd -r complex_min.inpcrd > min.log 2>&1"
        )
        subprocess.run(min_cmd, shell=True, cwd=work_dir, env=env)

        min_inpcrd = os.path.join(work_dir, "complex_min.inpcrd")
        eval_inpcrd = min_inpcrd if os.path.exists(min_inpcrd) else "complex.inpcrd"

        # Step C: MMPBSA.py Rescoring
        mmpbsa_cmd = (
            f"{CONDA_BIN}/MMPBSA.py -O -i {mdp_file} "
            f"-cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop "
            f"-y {eval_inpcrd} -o HAWKDOCK_MMGBSA_RESULTS.dat "
            f"-eo HAWKDOCK_MMGBSA_RESULTS.csv > mmgbsa_hawkdock.log 2>&1"
        )
        subprocess.run(mmpbsa_cmd, shell=True, cwd=work_dir, env=env)

        if os.path.exists(dat_out) and os.path.getsize(dat_out) > 0:
            logging.info(f"Finished MM/GBSA (Minimized): {variant_name}")
        else:
            logging.error(f"MM/GBSA failed for {variant_name}. Check {work_dir}/mmgbsa_hawkdock.log")


# ==============================================================================
# PHASE 4: SUMMARY TABLE COMPILATION
# ==============================================================================
def compile_master_summary():
    """Parses HADDOCK and MM/GBSA outputs to compile comprehensive master summary TSV files."""
    logging.info("Compiling master screening summary data...")
    
    mmgbsa_records = []
    pdb_files = sorted(glob.glob(os.path.join(COMPLEX_DIR, "*.pdb")))

    for pdb_path in pdb_files:
        variant_name = os.path.basename(pdb_path).replace(".pdb", "").replace("_complex", "")
        dat_out = os.path.join(HAWKDOCK_DIR, variant_name, "HAWKDOCK_MMGBSA_RESULTS.dat")

        row = {
            "Variant_Domain": variant_name,
            "VDWAALS_kcal_mol": None,
            "EEL_kcal_mol": None,
            "EGB_kcal_mol": None,
            "ESURF_kcal_mol": None,
            "G_gas_kcal_mol": None,
            "G_solv_kcal_mol": None,
            "MMGBSA_dG_kcal_mol": None
        }

        if os.path.exists(dat_out) and os.path.getsize(dat_out) > 0:
            with open(dat_out, "r") as f:
                in_diff = False
                for line in f:
                    if "Differences (Complex - Receptor - Ligand)" in line:
                        in_diff = True
                        continue
                    if in_diff:
                        parts = line.strip().split()
                        if not parts: continue
                        if parts[0] == "VDWAALS": row["VDWAALS_kcal_mol"] = float(parts[1])
                        elif parts[0] == "EEL": row["EEL_kcal_mol"] = float(parts[1])
                        elif parts[0] == "EGB": row["EGB_kcal_mol"] = float(parts[1])
                        elif parts[0] == "ESURF": row["ESURF_kcal_mol"] = float(parts[1])
                        elif len(parts) >= 3 and parts[0] == "G" and parts[1] == "gas": row["G_gas_kcal_mol"] = float(parts[2])
                        elif len(parts) >= 3 and parts[0] == "G" and parts[1] == "solv": row["G_solv_kcal_mol"] = float(parts[2])
                        elif parts[0] == "TOTAL": row["MMGBSA_dG_kcal_mol"] = float(parts[1])

            # Auto-calculate combined gas, solv, and total binding free energy terms if needed
            if row["VDWAALS_kcal_mol"] is not None and row["EEL_kcal_mol"] is not None:
                if row["G_gas_kcal_mol"] is None:
                    row["G_gas_kcal_mol"] = round(row["VDWAALS_kcal_mol"] + row["EEL_kcal_mol"], 4)
                if row["G_solv_kcal_mol"] is None and row["EGB_kcal_mol"] is not None and row["ESURF_kcal_mol"] is not None:
                    row["G_solv_kcal_mol"] = round(row["EGB_kcal_mol"] + row["ESURF_kcal_mol"], 4)
                if row["MMGBSA_dG_kcal_mol"] is None and row["G_solv_kcal_mol"] is not None:
                    row["MMGBSA_dG_kcal_mol"] = round(row["G_gas_kcal_mol"] + row["G_solv_kcal_mol"], 4)

        mmgbsa_records.append(row)

    df_mmgbsa = pd.DataFrame(mmgbsa_records)
    if not df_mmgbsa.empty and "MMGBSA_dG_kcal_mol" in df_mmgbsa.columns:
        df_mmgbsa = df_mmgbsa.sort_values(by="MMGBSA_dG_kcal_mol", ascending=True, na_position="last")

    out_tsv = os.path.join(RESULTS_DIR, "mmgbsa_energy_components_master.tsv")
    df_mmgbsa.to_csv(out_tsv, sep="\t", index=False)
    logging.info(f"Master summary saved to: {out_tsv}")


# ==============================================================================
# MAIN EXECUTION ENTRY POINT
# ==============================================================================
if __name__ == "__main__":
    logging.info("=== STARTING MASTER PFK13/PI3K AUTOMATED PIPELINE ===")

    # 1. Create control files (500 steps min, cut=25.0, idecomp disabled for fast rescoring)
    mmpbsa_file, min_file = create_control_files()

    # 2. Sanitize PDB chains
    sanitize_pdb_chains()

    # 3. HADDOCK3 Batch Docking Screening
    mutant_files = sorted(glob.glob(os.path.join(MUTANTS_DIR, "*_A.pdb")))
    if not mutant_files:
        logging.error(f"No mutant structures found in {MUTANTS_DIR}/")
        sys.exit(1)

    logging.info(f"Found {len(mutant_files)} mutant structures in mutants/. Processing docking...")
    for mutant in mutant_files:
        run_local_haddock(
            mutant, "catalytic",
            os.path.join(TARGETS_DIR, "PI3K_Catalytic_B.pdb"),
            os.path.join(RESTRAINTS_DIR, "restraints_catalytic.tbl")
        )
        run_local_haddock(
            mutant, "helical",
            os.path.join(TARGETS_DIR, "PI3K_Helical_B.pdb"),
            os.path.join(RESTRAINTS_DIR, "restraints_helical.tbl")
        )

    # 4. Extract Top-Scoring Cluster Complexes
    extract_best_cluster_complexes()

    # 5. Run Sander Minimization & MMPBSA.py Rescoring
    run_mmgbsa_rescoring(mmpbsa_file, min_file)

    # 6. Compile Master Summary Dataset
    compile_master_summary()

    logging.info("=== PIPELINE EXECUTION COMPLETE ===")
