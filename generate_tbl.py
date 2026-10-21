import os

# Active residues mapping
active_pfk13 = [135, 137, 184, 230, 232, 278, 327, 372]
active_pi3k_cat = [30, 40, 66, 80, 145, 158]
active_pi3k_hel = [37, 59, 61, 63, 93, 94, 96, 97]

def write_tbl(filename, active_b):
    with open(filename, "w") as f:
        # Force PfK13 active sites (Segid A) to contact PI3K active sites (Segid B)
        for res_a in active_pfk13:
            f.write(f"assign (resi {res_a} and segid A)\n       (\n")
            b_lines = [f"        (resi {res_b} and segid B)" for res_b in active_b]
            f.write("\n     or ".join(b_lines) + "\n       ) 2.0 2.0 0.0\n\n")

        # Force PI3K active sites to contact PfK13 active sites
        for res_b in active_b:
            f.write(f"assign (resi {res_b} and segid B)\n       (\n")
            a_lines = [f"        (resi {res_a} and segid A)" for res_a in active_pfk13]
            f.write("\n     or ".join(a_lines) + "\n       ) 2.0 2.0 0.0\n\n")

os.makedirs("restraints", exist_ok=True)
write_tbl("restraints/restraints_catalytic.tbl", active_pi3k_cat)
write_tbl("restraints/restraints_helical.tbl", active_pi3k_hel)
