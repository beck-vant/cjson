import MDAnalysis as mda
from ase.io import read
from ase import Atoms
import numpy as np
import json
import cjson


# H2 = Atoms(["H"] * 2)
# print(H2)


# print(mol)

# mdamol = cjson.chemical_json_to_mda(mol)

# print(mdamol)

# mol = cjson.atomgroup_to_cjson(
#     mdamol.select_atoms("all"),
# )

U = mda.Universe("peptide.pdb")
print("Original:", U)
print(U.select_atoms("all").chainIDs)
mol = cjson.atomgroup_to_cjson(
    U.select_atoms("all"),
)

with open("peptide.cjson", "w") as jout:
    jout.write(cjson.dump(mol))
    # json.dump(mol, jout, sort_keys=False, indent=2)

mdamol = cjson.chemical_json_to_mda(mol)
mdamol.select_atoms("all").write("Recreated.pdb")
print("Recreated:", mdamol)


asemol = cjson.chemical_json_to_ase(mol)
asemol.write("Recreated_ase.pdb")
