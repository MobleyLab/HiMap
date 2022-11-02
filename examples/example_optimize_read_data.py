import himap

"""
This example file walks through how to read in scores such as similarity and
IDs from external files to then determine optimal FEP designs.
This uses an interactive clustering mode. To select defaults during
clustering, hit enter or type the displayed options.
The user selected clusters will then be sent for design optimization.
"""
# *****************************************************************************
# This example file was written by Dr. Mary Pitman. 2022
# *****************************************************************************

#-------------------------------------------------------#
# Define input files, read data.
#-------------------------------------------------------#
# Input files for weight scores and ligand names.
sim_scores_in = '../test/optimize/sim_scores.csv'
IDs_in = '../test/optimize/mol_names.txt'

# Read files, clean any potential NaN scores.
#   Added optional parameter:
#             delimiter: default is ','
n_arr, ID_list = himap.read_data(sim_scores_in, IDs = IDs_in)

#-------------------------------------------------------#
# Clustering.
#-------------------------------------------------------#
# Perform clustering.
#   sub_arr, sub_ID:   the n_arr and ID_list subdivided by clusters
#   selected_clusters: user selected clusters during interaction.
sub_arr, sub_ID, selected_clusters = himap.cluster_interactive(n_arr, ID_list)

#-------------------------------------------------------#
# Optimization.
#-------------------------------------------------------#
# Example reference ligands.
ref_ligs = ['mol_0', 'mol_1', 'mol_2', 'mol_3', 'mol_4']

# Send the user selected clusters for optimization.
# Send the user selected clusters for optimization.
himap.clusters2optimize(sub_arr, sub_ID, clusters2optim = selected_clusters,
                        ref_ligs=ref_ligs, num_edges = '2n', optim_types = ['A', 'D']
                       )
