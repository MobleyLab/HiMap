## Contents of this directory

 * [clustering.py](https://github.com/MobleyLab/HiMap/tree/main/himap#clustering)
 * [optimal.py](https://github.com/MobleyLab/HiMap/tree/main/himap#initiate-optimization-in-python)
 * [optimal_design.R](https://github.com/MobleyLab/HiMap/tree/main/himap#run-optimization)
 * [utils.py](https://github.com/MobleyLab/HiMap/tree/main/himap#utilities)


Clustering
-------

Functions for if user wants to cluster ligands and
    assign the clustering neighbor distance cutoff. Uses DBSCAN.
    
        Parameters:
            sim_data: the similarity array from similarity scores
            ID_list: list of ID names
             
        Returns:
            sub_arrs: dict, n_arr subdivided into dict of cluster number keys
            sub_IDs: dict, ID_list subdivided into dict of cluster number keys
            selected_clusters: controls optimization feature for which clusters
                           are run in optimization.
                        Options are 'all',
                        Run only clusters with reference ligands, 'w_ref_lig'
                        Run only clusters numerically specified by user_clusters,
                        in other words, user interactive input or numeric arguments.
    

Minimal example of usage: 

```python
import himap

# Perform clustering.
#   sub_arr, sub_ID:   the n_arr and ID_list subdivided by clusters
#   selected_clusters: user selected clusters during interaction.
sub_arr, sub_ID, selected_clusters = himap.cluster_interactive(n_arr, ID_list)

# Send the user selected clusters for optimization.
himap.clusters2optimize(sub_arr, sub_ID, clusters2optim = selected_clusters)
```

Initiate Optimization in Python
-------

Either initiates optimization without running clustering or
processed clustered data to run optimization. 
    
        Parameters:
            np_arr = 2D similarity array, for example
                strict_numpy
                    
        Optional Parameters:
            optim_types = [ optim1 , optim2 ](str), optimization type such as
                        'A' and 'D'. Currently two types are required. Options are:
                        'A', 'D', 'P', 'mA', 'mP', 'negA', 'negD', 'random'
            
            db_mol = output of lomap.DBMolecules()
                
            ref_lig (list of str or str) = user input reference ligand.
                If not input, will be calculated based on max similarity.
                
            ID_list (str) = the user can define ligand names. The default
                is not not enter a list but instead output.
                
            num_edges = the number of edges requested for optimization. If 'n' is the
                        number of ligands in each cluster, the options are:
                        '1n', '2n', 'nlnn', 'min', 'max', integers.
                        The edge number requested must be in the range of [min, max].
                        If less than min, will be set to min. If greater than max,
                        will be set to max.
        
        Returns:
            Optimal graph outputs.
                    
        Example usage:
        himap.Optimize(strict_numpy, db_mol, ref_lig = 'lig_1')

Run Optimization
-----
Performs graph optimization using the Fedorov Exchange algorithm. 
This file:

  * runs the optimization via Fedorov.exchange
  * computes the optimal designs with and without weighting,
  * and outputs the optimal designs.

Utilities
-----
Utility functions for handling clustering and optimization. Includes: 
  * data reading
  * generation of random similarity scores
  * cleaning up invalid data values or ligand names
  * writing CSV files for export
  * flexible processing of strings for reading interactive user inputs. 
