import json
import numpy as np
from MDAnalysis.core.groups import AtomGroup
import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import (
    Atomnames,
    Elements,
    Resids,
    Resnames,
    Segids,
    Masses,
    Charges,
)
from ase.data import atomic_numbers, chemical_symbols

symbol_to_number = atomic_numbers  # {'H': 1, 'He': 2, ...}
number_to_symbol = chemical_symbols  # ['X', 'H', 'He', ...]


def chemical_json_to_mda(cjson_data: dict) -> mda.Universe:
    """Create MDAnalysis Universe from Chemical JSON with proper chainID handling.

    Args:
        cjson_data: Dictionary in Chemical JSON format

    Returns:
        MDAnalysis Universe object
    """
    atoms_data = cjson_data["atoms"]
    elements = np.array(atoms_data["elements"]["number"], dtype=np.int32)
    elements = [number_to_symbol[x] for x in elements]
    n_atoms = len(elements)

    # Get chainIDs if available
    chainIDs = atoms_data.get("chainIDs", ["A"] * n_atoms)

    # Process residue information
    residue_data = cjson_data.get("properties", {}).get("residues", [])
    n_residues = len(residue_data) if residue_data else 1

    # Prepare topology arrays
    resindex = np.zeros(n_atoms, dtype=np.int32)
    resids = []
    resnames = []
    segids = []

    for i, res in enumerate(residue_data):
        for atom_idx in res["atoms"]:
            resindex[atom_idx] = i
        resids.append(res.get("resid", i + 1))
        resnames.append(res.get("resname", "UNK"))
        segids.append(res.get("segid", "SYST"))

    # Create segments based on chainIDs
    unique_segments = list(set(chainIDs)) or ["A"]

    # Create mapping from residue index to segment index
    residue_segindex = []
    if residue_data:
        # For each residue, get the chainID of its first atom
        for res in residue_data:
            if res["atoms"]:  # If residue has atoms
                first_atom_idx = res["atoms"][0]
                chainID = chainIDs[first_atom_idx]
                residue_segindex.append(unique_segments.index(chainID))
            else:
                residue_segindex.append(0)  # Default to first segment
    else:
        # If no residue data, assign all atoms to first segment
        residue_segindex = [0] * n_residues

    # Create universe
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_residues,
        n_segments=len(unique_segments),
        atom_resindex=resindex,
        residue_segindex=residue_segindex,  # Now correctly sized
        trajectory=True,
    )

    # Set basic attributes
    u.add_TopologyAttr(
        "name", atoms_data.get("labels", [str(i) for i in range(n_atoms)])
    )
    u.add_TopologyAttr("type", ["X"] * n_atoms)
    u.add_TopologyAttr("element", elements)
    u.atoms.positions = np.array(atoms_data["coords"]["3d"]).reshape(-1, 3)

    # Set masses/charges
    if "masses" in atoms_data:
        u.add_TopologyAttr("mass", np.array(atoms_data["masses"], dtype=np.float32))
    if "formalCharges" in atoms_data:
        u.add_TopologyAttr(
            "charge", np.array(atoms_data["formalCharges"], dtype=np.float32)
        )

    # Set residue and segment attributes
    u.add_TopologyAttr("resid", np.array(resids, dtype=np.int32))
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("segid", unique_segments)

    # Add chainIDs as a topology attribute
    u.add_TopologyAttr("chainIDs", chainIDs)

    # Add bonds if present
    if "bonds" in cjson_data:
        connections = np.array(
            cjson_data["bonds"]["connections"]["index"], dtype=np.int32
        )
        u.add_TopologyAttr("bonds", connections.reshape(-1, 2))

    # Add unit cell if present
    if "unitCell" in cjson_data:
        cell = cjson_data["unitCell"]
        if "cellVectors" in cell:
            u.dimensions = np.array(cell["cellVectors"], dtype=np.float32).reshape(3, 3)
        else:
            u.dimensions = np.array(
                [
                    cell["a"],
                    cell["b"],
                    cell["c"],
                    cell["alpha"],
                    cell["beta"],
                    cell["gamma"],
                ],
                dtype=np.float32,
            )

    return u


def atomgroup_to_cjson(ag: AtomGroup, name: str = None) -> dict:
    """Convert MDAnalysis AtomGroup to Chemical JSON format.

    Args:
        ag: MDAnalysis AtomGroup (selection)
        name: Optional name for the system

    Returns:
        Dictionary compliant with Chemical JSON schema
    """

    # Basic structure
    cjson = {
        "chemicalJson": 1,
        "atoms": {
            "elements": {"number": [symbol_to_number[x] for x in ag.elements]},
            "coords": {
                "3d": ag.positions.flatten().tolist()  # Flatten to [x1,y1,z1, x2,y2,z2,...]
            },
            "labels": ag.names.tolist(),  # Atom names (e.g., "CA", "N")
            "chain": (
                ag.chainIDs.tolist() if hasattr(ag, "chainIDs") else ["A"] * len(ag)
            ),
        },
        "name": name if name else "MDAnalysis Selection",
        # "formula": ag.universe.atoms[:].guess_chemical_formula(),
    }

    # Add residue information as properties
    residue_data = []
    for res in ag.residues:
        entry = {
            "resid": int(res.resid) if hasattr(res, "resid") else 1,
            "resname": str(res.resname) if hasattr(res, "resname") else "UNK",
            "segid": str(res.segid) if hasattr(res, "segid") else "SYST",
            "atoms": [int(i) for i in res.atoms.indices],
        }

        residue_data.append(entry)

    cjson["properties"] = {
        "residues": residue_data,
        "totalCharge": 0,  # Default, modify if needed
        "nAtoms": ag.n_atoms,
        "nResidues": ag.n_residues,
    }

    # Add bonds if available (e.g., from PSF/TOPology)
    if hasattr(ag, "bonds"):
        bond_connections = []
        bond_orders = []
        for bond in ag.bonds:
            bond_connections.extend([bond.atoms[0].ix, bond.atoms[1].ix])
            bond_orders.append(1)  # Default order=1 (single bond)

        cjson["bonds"] = {
            "connections": {"index": bond_connections},
            "order": bond_orders,
        }

    # Add unit cell if periodic
    if ag.universe.dimensions is not None:
        a, b, c, alpha, beta, gamma = ag.universe.dimensions[:6]
        cjson["unitCell"] = {
            "a": a,
            "b": b,
            "c": c,
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "cellVectors": ag.universe.trajectory.ts.triclinic_dimensions.flatten().tolist(),
        }

    return cjson
