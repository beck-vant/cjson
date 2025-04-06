from ase import Atoms
import numpy as np
import json


def ase_to_chemical_json(atoms: Atoms, name: str = None) -> dict:
    """Convert ASE Atoms object to Chemical JSON format.

    Args:
        atoms: ASE Atoms object
        name: Optional name for the system

    Returns:
        Dictionary compliant with Chemical JSON schema
    """
    # Basic structure
    cjson = {
        "chemicalJson": 1,
        "name": name if name else "ASE System",
        "formula": atoms.get_chemical_formula(),
        "atoms": {
            "elements": {"number": atoms.get_atomic_numbers().tolist()},
            "coords": {
                "3d": atoms.get_positions()
                .flatten()
                .tolist()  # Flatten to [x1,y1,z1, x2,y2,z2,...]
            },
        },
        "properties": {
            "totalCharge": (
                int(atoms.get_initial_charges().sum())
                if atoms.has("initial_charges")
                else 0
            ),
            "molecularMass": float(atoms.get_masses().sum()),
        },
    }

    # Add cell information if periodic
    if any(atoms.pbc):
        cell = atoms.get_cell()
        cjson["unitCell"] = {
            "a": float(np.linalg.norm(cell[0])),
            "b": float(np.linalg.norm(cell[1])),
            "c": float(np.linalg.norm(cell[2])),
            "alpha": float(
                np.degrees(
                    np.arccos(
                        np.dot(cell[1], cell[2])
                        / (np.linalg.norm(cell[1]) * np.linalg.norm(cell[2]))
                    )
                )
            ),
            "beta": float(
                np.degrees(
                    np.arccos(
                        np.dot(cell[0], cell[2])
                        / (np.linalg.norm(cell[0]) * np.linalg.norm(cell[2]))
                    )
                )
            ),
            "gamma": float(
                np.degrees(
                    np.arccos(
                        np.dot(cell[0], cell[1])
                        / (np.linalg.norm(cell[0]) * np.linalg.norm(cell[1]))
                    )
                )
            ),
            "cellVectors": cell.flatten().tolist(),
        }

    # Add velocities if present
    if atoms.has("velocities"):
        cjson["atoms"]["coords"]["3dSets"] = [
            atoms.get_positions().tolist(),
            atoms.get_velocities().tolist(),
        ]

    # Add momenta if present
    if atoms.has("momenta"):
        cjson["properties"]["momenta"] = atoms.get_momenta().flatten().tolist()

    # Add charges if present
    if atoms.has("initial_charges"):
        cjson["atoms"]["formalCharges"] = atoms.get_initial_charges().tolist()

    # Add calculator results if present
    if atoms.calc is not None and hasattr(atoms.calc, "results"):
        calc_results = atoms.calc.results
        if "energy" in calc_results:
            cjson["properties"]["totalEnergy"] = float(calc_results["energy"])
        if "forces" in calc_results:
            cjson["properties"]["forces"] = calc_results["forces"].flatten().tolist()
        if "stress" in calc_results:
            cjson["properties"]["stress"] = calc_results["stress"].flatten().tolist()
        if "dipole" in calc_results:
            cjson["properties"]["dipoleMoment"] = calc_results["dipole"].tolist()

    return cjson


def chemical_json_to_ase(cjson_data: dict) -> Atoms:
    """Create an ASE Atoms object from Chemical JSON data.

    Args:
        cjson_data: Dictionary in Chemical JSON format

    Returns:
        ASE Atoms object with all available information
    """
    assert type(cjson_data) is dict, "pass the object not the string"
    atoms_data = cjson_data["atoms"]

    # Extract basic atom information
    atomic_numbers = np.array(atoms_data["elements"]["number"], dtype=int)
    positions = np.array(atoms_data["coords"]["3d"]).reshape(-1, 3)

    # Create basic Atoms object
    atoms = Atoms(
        numbers=atomic_numbers,
        positions=positions,
        pbc=False,  # Will update if unit cell exists
    )

    # Set atom names if available
    if "labels" in atoms_data:
        atoms.set_array("names", np.array(atoms_data["labels"]))

    # Set masses if available
    if "masses" in atoms_data:
        atoms.set_masses(np.array(atoms_data["masses"], dtype=float))

    # Set charges if available
    if "formalCharges" in atoms_data:
        atoms.set_charges(np.array(atoms_data["formalCharges"], dtype=float))

    # Set velocities if available (in 3dSets)
    if (
        "coords" in atoms_data
        and "3dSets" in atoms_data["coords"]
        and len(atoms_data["coords"]["3dSets"]) > 1
    ):
        velocities = np.array(atoms_data["coords"]["3dSets"][1])
        atoms.set_velocities(velocities)

    # Set residue information if available
    if "properties" in cjson_data and "residues" in cjson_data["properties"]:
        residue_info = []
        for res in cjson_data["properties"]["residues"]:
            for atom_idx in res["atoms"]:
                residue_info.append(
                    {
                        "resid": res.get("resid", 1),
                        "resname": res.get("resname", "UNK"),
                        "chainID": res.get("chainID", "A"),
                    }
                )
        # Convert to numpy arrays
        resids = np.array([r["resid"] for r in residue_info])
        resnames = np.array([r["resname"] for r in residue_info])
        chainIDs = np.array([r["chainID"] for r in residue_info])

        atoms.set_array("resid", resids)
        atoms.set_array("resname", resnames)
        atoms.set_array("chainID", chainIDs)

    # Set unit cell if available
    if "unitCell" in cjson_data:
        cell = cjson_data["unitCell"]
        if "cellVectors" in cell:
            cell_array = np.array(cell["cellVectors"], dtype=float).reshape(3, 3)
            atoms.set_cell(cell_array)
            atoms.set_pbc(True)
        else:
            cell_params = [
                cell["a"],
                cell["b"],
                cell["c"],
                cell["alpha"],
                cell["beta"],
                cell["gamma"],
            ]
            atoms.set_cell(cell_params)
            atoms.set_pbc(True)

    # Set calculator results if available
    if "properties" in cjson_data:
        props = cjson_data["properties"]
        if "totalEnergy" in props:
            atoms.calc.results["energy"] = float(props["totalEnergy"])
        if "forces" in props:
            atoms.calc.results["forces"] = np.array(props["forces"]).reshape(-1, 3)
        if "stress" in props:
            atoms.calc.results["stress"] = np.array(props["stress"])

    return atoms
