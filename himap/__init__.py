"""
HiMap
======

Alchemical free energy calculations hold increasing promise as an aid to drug 
discovery efforts. However, applications of these techniques in discovery 
projects have been relatively few, partly because of the difficulty of planning 
and setting up calculations. The High Informaiton Mapper (HiMap) is an
automated algorithm to optimize relative free energy calculations between
potential ligands within a substantial of compounds.

Authors: Mary Pitman      <mpitman@uci.edu>
         David Mobley     <dmobley@uci.edu>
         
Licence: GNU 3.0

Using
-----
import lomap
import himap

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
n_arr = himap.clean_NaN(sim_np)

#-------------------------------------------------------#
# Clustering.
#-------------------------------------------------------#
# Create ID_list from db_mol prior to clustering.
ID_list = lomap.db_mol_IDs(db_mol, n_arr)

# Perform clustering.
#   sub_arr, sub_ID:   the n_arr and ID_list subdivided by clusters
#   selected_clusters: user selected clusters during interaction.
sub_arr, sub_ID, selected_clusters = himap.cluster_interactive(n_arr, ID_list)

#-------------------------------------------------------#
# Optimization.
#-------------------------------------------------------#
# Example reference ligand.
ref_ligs = ['ejm_31']

# Send the user selected clusters for optimization.
himap.clusters2optimize(sub_arr, sub_ID, clusters2optim = selected_clusters, ref_ligs=ref_ligs)
"""

from .optimal import Optimize
from .utils import *
from .clustering import *

