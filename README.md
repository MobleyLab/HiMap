[![CI](https://github.com/OpenFreeEnergy/Lomap/actions/workflows/CI.yaml/badge.svg)](https://github.com/OpenFreeEnergy/Lomap/actions/workflows/CI.yaml)
## Contents of this file

 * HiMap Introduction
 * Requirements
 * Authors
 * Installation
 * Usage
 * Troubleshooting

HiMap Introduction
-------

HiMap includes design generation based on statistical optimality. 
Alchemical free energy calculations hold increasing promise 
as an aid to drug discovery efforts. However, applications of 
these techniques in discovery projects have been relatively 
rare, partly because of the difficulty of planning and setting up 
calculations. The lead optimization mapper (LOMAP) was 
introduced as an automated algorithm to plan relative 
free energy calculations between potential ligands. LOMAP was further
developed to be based on free, open-source APIs such as RDKit. HiMap
now includes clustering of ligand series, and optimization of free
energy perturbation networks. 

Requirements
-------
* RDKit Release > 2021
* NetworkX
* Matplotlib 
* python >= 3.8
* R
* rpy2=3.4.5
* kneed=0.7.0
* scikit-learn=0.23.2
* scipy


Authors
-------

Contact for the graph optimization features within this fork:
* Mary Pitman <mpitman@uci.edu>
* David Mobley <dmobley@uci.edu>
    
Predecessor authors:  
* Gaetano Calabro' <gcalabro@uci.edu>
* Mark Mackey
* Lester Hedges
* Antonia S J S Mey
* Jenke Scheen

Installation
-----

For optimal lomap methods, build the conda environment and install Lomap from file:

https://github.com/pitmanme/HiMap/blob/main/devtools/conda-envs/himap_env.yml

with:

`conda env create -f himap_env.yml python=3.8`


Usage
-----
#### Example scripts are included for various purposes:
* To run HiMap with optimization \
    `python examples/example_optimize.py`
* To read in scores and optimize \
    `python examples/example_optimize_read_data.py`
* For a basic LOMAP run without optimization \
    `python examples/example.py`
* For generating radial graphs with a hub \
    `python examples/example_radial.py`

#### To run LOMAP without optimization, as a commandline tool:
`
lomap test/basic/
`


#### If you would rather use the API directly:
* To generate optimal designs, try:

```python
import lomap

#-------------------------------------------------------#
# Generate similarity scores.
#-------------------------------------------------------#
# Read molecules from test directory.
db_mol = lomap.DBMolecules('../test/radial/', output=True, radial=True)
    
# Generate the strict and loose symmetric similarity score matrices.
strict, loose = db_mol.build_matrices()
    
# Convert the similarity matrix to numpy array
sim_np = strict.to_numpy_2D_array()

# Clean data if Lomap produces rare error. If score is NaN, replace with 0.0
n_arr = lomap.clean_NaN(sim_np)

#-------------------------------------------------------#
# Clustering.
#-------------------------------------------------------#
# Create ID_list from db_mol prior to clustering.
ID_list = lomap.db_mol_IDs(db_mol, n_arr)

# Perform clustering.
#   sub_arr, sub_ID:   the n_arr and ID_list subdivided by clusters
#   selected_clusters: user selected clusters during interaction.
sub_arr, sub_ID, selected_clusters = lomap.cluster_interactive(n_arr, ID_list)

#-------------------------------------------------------#
# Optimization.
#-------------------------------------------------------#
# Example reference ligand.
ref_ligs = ['ejm_31']

# Send the user selected clusters for optimization.
lomap.clusters2optimize(sub_arr, sub_ID, clusters2optim = selected_clusters, ref_ligs=ref_ligs)
```


* To generate optimal designs using external scores or weights, try:

```python
import lomap

#-------------------------------------------------------#
# Define input files, read data.
#-------------------------------------------------------#
# Input files for weight scores and ligand names.
sim_scores_in = '../test/optimize/sim_scores.csv'
IDs_in = '../test/optimize/mol_names.txt'

# Read files, clean any potential NaN scores.
#   Added optional parameter:
#             delimiter: default is ','
n_arr, ID_list = lomap.read_data(sim_scores_in, IDs = IDs_in)

#-------------------------------------------------------#
# Clustering.
#-------------------------------------------------------#
# Perform clustering.
#   sub_arr, sub_ID:   the n_arr and ID_list subdivided by clusters
#   selected_clusters: user selected clusters during interaction.
sub_arr, sub_ID, selected_clusters = lomap.cluster_interactive(n_arr, ID_list)

#-------------------------------------------------------#
# Optimization.
#-------------------------------------------------------#
# Example reference ligands.
ref_ligs = ['mol_0', 'mol_1', 'mol_2', 'mol_3', 'mol_4']

# Send the user selected clusters for optimization.
lomap.clusters2optimize(sub_arr, sub_ID, clusters2optim = selected_clusters,
                        ref_ligs=ref_ligs, num_edges = '2n', optim_types = ['A', 'D']
                        )
```


* To generate original LOMAP designs, try:

```python
import lomap

# Generate the molecule database starting from a directory containing .mol2 files

db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True)

    #More graphing options:
    # Use the complete radial graph option. The ligand with the most structural similarity to all of the others will be picked as the 'lead compounds' and used as the central compound.
    db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, radial=True)

    # Use a radial graph with a manually specified hub compound
    db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, radial=True, hub=filename.mol2)

    # Use a radial graph with a manually specified hub compound and fast graphing option
    #the fast graphing option create the initial graph by connecting the hub ligand with the possible surrounding ligands and add surrounding edges based on the similarities accoss surrounding nodes
    db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, radial=True, hub=filename.mol2, fast=True)

# Calculate the similarity matrix betweeen the database molecules. Two molecules are generated
# related to the scrict rule and loose rule 

strict, loose = db_mol.build_matrices()

# Generate the NetworkX graph and output the results
nx_graph = db_mol.build_graph() 


# Calculate the Maximum Common Subgraph (MCS) between 
# the first two molecules in the molecule database 
# ignoring hydrogens and depicting the mapping in a file
    
MC = lomap.MCS.getMapping(db_mol[0].getMolecule(), db_mol[1].getMolecule(), hydrogens=False, fname='mcs.png')


# Alchemical transformation are usually performed between molecules with
# the same charges. However, it is possible to allow this transformation
# manually setting the electrostatic score for the whole set of molecules 
# producing a connected graph. The electrostatic scrore must be in the 
# range [0,1]


db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, ecrscore=0.1)
strict, loose = db_mol.build_matrices()
nx_graph = db_mol.build_graph() 
```

Troubleshooting
-----
* Why is optimization not finding a random seed design? \
Check that similarity scores are not zero, nearly zero, or in large part near zero. If the similarity scores are non-zero but very dissimilar, try decreasing the neighbor distance cutoff value (epsilon) selected during clustering. If the similarity scores are zero or near zero, this may be an example of a rare limitation, currently in development, in the LOMAP similarity scores. Check that the ligands are indeed very dissimilar. If they are not, consider using an alternate similarity metric for the time being.
