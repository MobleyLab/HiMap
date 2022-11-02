# ******************
# MODULE DOCSTRING
# ******************

"""

LOMAP: Maximum Common Subgraph and scoring calculations
=====

Alchemical free energy calculations hold increasing promise as an aid to drug
discovery efforts. However, applications of these techniques in discovery
projects have been relatively few, partly because of the difficulty of planning
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an
automated algorithm to plan efficient relative free energy calculations between
potential ligands within a substantial of compounds. The optimal module
generates A and D optimal graphs to reduce uncertainty in FE results for
RBFE calculations.

"""

# *****************************************************************************
# This module was written by Dr. Mary Pitman. 2022
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see http://www.gnu.org/licenses/
# *****************************************************************************


# ****************
# MODULE IMPORTS
# ****************
import numpy as np
import pandas as pd
import os
import lomap

# Imports specific to python to R interaction:
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import globalenv

# Must activate rpy2.
pandas2ri.activate()

__all__ = ['df_gen', 'ref_lig_gen', 'Optimize']

# ****************
# FUNCTIONS
# ****************

def df_gen(np_arr, **kwargs):
    ''' Pandas dataframe generation from lomap
    generated object db_mol and numpy 2D symmetric
    array of chemical similarities.
    
        Parameters:
            np_arr : 2D similarity array
            db_mol : output of lomap.DBMolecules()
        
       Returns:
            df : pandas dataframe of similarity scores
                 with ligands names as row and col IDs.
    '''
    # Read in option input variables.
    ID_list = kwargs.get('ID_list', None)
    db_mol = kwargs.get('db_mol', None)
    
    # If no list of IDs, create one.
    if db_mol is None:
        if ID_list is None:
            mol_names = ['ID']
            for i in range(np_arr.shape[1]):
                fname = i
                mol_names.append(fname)
            # Convert the name list to an np array.
            mol_arr_pre = np.array(list(mol_names))
            # Cleave the 'ID' entry
            mol_arr = mol_arr_pre[1:]
        else:
            mol_arr = np.array(list(ID_list))
    # If similarity was not read in.
    if db_mol is not None:
        if ID_list is None:
            # Generate a list of the ligand names.
            mol_names = ['ID']
            for i in range(np_arr.shape[1]):
                fname = db_mol[i].getName()
                fsplit = fname.split(".", 1)[0]
                mol_names.append(fsplit)
            # Convert the name list to an np array.
            mol_arr_pre = np.array(list(mol_names))
            # Cleave the 'ID' entry
            mol_arr = mol_arr_pre[1:]
        else:
            mol_arr = np.array(list(ID_list))
        
    # Create the pandas DataFrame to pass to R.
    df = pd.DataFrame(np_arr)

    # Set the df titles to the ligand IDs
    # First entry "ID" is ommitted.
    df.set_axis(mol_arr, axis=0, inplace = True)
    df.set_axis(mol_arr, axis=1, inplace = True)
    return df
    

def py_run_optimization(*args):
    '''
    Runs function 'run_optimization' in lomap/optimal_design.R.
    Sends the R output back to python.
        
        Parameters:
            ref_lig: reference ligand to use, selected or calculated
            dataframe : r dataframe generated from similarity scores
            optim_type1 and optim_type2: optimization type such as
                        'A' and 'D'. Currently two types are required.
                        
        Returns:
            Optimal graph outputs
    '''
    # Get the directory of this script
    optim_dir = os.path.dirname(os.path.realpath(__file__))
    # Does this source other optimal_design.R if I am in a dif folder? 
    r = robjects.r
    # Need to find the R script and put it in. Don't know why this happened
    r.source('{}/optimal_design.R'.format(optim_dir))
    #r.source('optimal_design_backup.R')
    if args:
        out=r.run_optimization(*args)
    else:
        out=r.run_optimization()
    return out


def ref_lig_gen(df):
    ''' This function selects a reference ligand.
    Parameters:
    
        df: 2D symmetric pandas dataframe with similarity scores.
            Returned from optimal.df_gen().
            
    Returns:
        ref_lig: The reference ligand for optimal graph generation.
                 The ligand with the maximal sum similarity score.
    '''
    # Sum each row of the ligand similarity dataframe, df.
    df_summed = pd.DataFrame(np.sum(df.values, axis=1), columns=['sum'])
    # Find the index of the maximum value in summed df.
    x = df_summed.idxmax()
    # Retreive from original dataframe the ligand name at max sum.
    ref_lig = df.index[x[0]]
    return ref_lig


def Optimize(np_arr, **kwargs):
    ''' Main function.
    
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
        lomap.Optimize(strict_numpy, db_mol, ref_lig = 'lig_1')
    '''
    # Set and get optional input arguments. Defaults: lomap,
    # A/D Optimization, no reference ligand given.
    db_mol = kwargs.get('db_mol', None)
    optim_types = kwargs.get('optim_types', ['A', 'D'])
    ref_lig = kwargs.get('ref_lig', None)
    ID_list = kwargs.get('ID_list', None)
    num_edges = kwargs.get('num_edges', 'nlnn')
    
    # db_mol will be None if sim not calculated in LOMAP.
    if db_mol is None:
        # If a list of ligand names were not read in.
        if ID_list is None:
            df = df_gen(np_arr)
        else:
            df = df_gen(np_arr, ID_list = ID_list)
    if db_mol is not None:
        if ID_list is None:
            df = df_gen(np_arr, db_mol = db_mol)
        else:
            df = df_gen(np_arr, db_mol = db_mol, ID_list = ID_list)
    
    # Define the options for number of edges.
    # Get length of ligands to be optimized
    n = np_arr.shape[1]
    min_connect = n-1
    # Subtract 1 because ref ligs will be added.
    max_connect = (n*(n - 1)//2) - 1
    # Calculate edge number selection
    try:
        int(num_edges) == num_edges
        num_edges = num_edges
    except:
        if num_edges is 'nlnn':
            num_edges = round(n*np.log(n))
        elif num_edges is '1n':
            num_edges = n
        elif num_edges is '2n':
            num_edges = 2*n
        elif num_edges is 'min':
            num_edges = min_connect
        elif num_edges is 'max':
            num_edges = max_connect
        else:
            raise ValueError(f"Invalid num_edges input. Valid inputs are: "
                              "'nlnn', '1n', '2n', 'min', 'max', and int")
            
    # Test edge selection to ensure in range.
    # Design must be fully connected.
    if num_edges < min_connect or num_edges > max_connect:
        print(f"Requested edge number, {num_edges}, is out of bounds. "
               "Range is [{min_connect}, {max_connect}]")
        # If less than range make min_connect
        if num_edges < min_connect:
            num_edges = min_connect
            print(f"Preparing optimization with {min_connect} edges")
        # if more than range make max_connect
        if num_edges > max_connect:
            num_edges = max_connect
            print(f"Preparing optimization with {max_connect} edges")
    
    # For testing
    print(f"The input number of edges is {num_edges}")
    
    # Convert pandas df to R df.
    r_df = pandas2ri.py2rpy(df)
    # Convert optimization types list into R str vector
    r_optim_types = robjects.vectors.StrVector(optim_types)
    # Select the reference ligand for graph generation if none given.
    if ref_lig is None:
        ref_lig = ref_lig_gen(df)
    print("The reference ligand is", ref_lig)
    # Ouput optimal graphs using optimal_design.R.
    c=py_run_optimization(ref_lig, r_df, r_optim_types, num_edges)
