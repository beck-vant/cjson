# cjson
Chemical JSON tools for MM

Integration for [OpenChemistry Chemical JSON](https://github.com/OpenChemistry/chemicaljson/tree/main) to various python modules.

Currently supports:

Currently supports:

| Library      | Read | Write |
|--------------|------|-------|
| MDAnalysis   | ✔️   | ✔️    |
| ASE          | ✔️   | ✔️    |
| rdkit        | ❌   | ❌    |
| mdtraj       | ❌   | ❌    |


Chemical JSON is for a chemical, not a trajectory - it contains a single set of coordinates.

chain ID is NOT part of the original CJSON but we are adding it.