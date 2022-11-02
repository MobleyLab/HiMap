# ******************
# MODULE DOCSTRING
# ******************

"""

HiMap
=====

Alchemical free energy calculations hold increasing promise as an aid to drug
discovery efforts. However, applications of these techniques in discovery
projects have been relatively few, partly because of the difficulty of planning
and setting up calculations. The High Information Mapper (HiMap) is an
automated algorithm to optimize efficient relative free energy calculations between
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
import csv
import numpy as np
import sys
import json

__all__ = ['read_data', 'rand_sim_scores', 'clean_up', 'record_dicts', 'write_csv', 'clean_NaN', 'db_mol_IDs', 'multi_delim']

# ****************
# FUNCTIONS
# ****************

def clean_NaN(arr):
    '''
    If entry is NaN, replace with 0.0.
        
        Parameters:
            arr: array
        Returns:
            arr_cleaned: array with NaN replaced with 0.0
    '''
    arr_cleaned = np.where(np.isnan(arr), 0.0, arr)
    return arr_cleaned
    
    
def multi_delim(input):
    '''
    Cleans multiple delimiters to avoid read issues in
    csv or user inputs.
        Parameters:
            input: string that could have varied delimiters

        Returns:
            cleaned: input where only delimiter is " ".
    '''
    delimiters = [" ", ",", "[", "]", ")", "("]
    cleaned = input
    for i in delimiters:
        cleaned = cleaned.replace(i, " ")
    return cleaned


def read_data(sim_scores, **kwargs):
    '''
    Reads in similarity scores and ligand IDs from
    an external file if called.
    
        Parameters:
            sim_scores: a CSV file
            
        Optional Parameters:
            IDs: a txt or plain text file,
                Default: None, if no IDs are given,
                will generate number list
            delimiter: the default is ',' you can
                provide other delimiters as needed.
        
        Returns:
            arr_cleaned: numpy array of similarity scores
            ID_list: list of ID names
    
    Example usage: read_data(sim_scores, IDs = IDs)
    '''
    # Change to make: allow various delimiters
    
    # Define optional arugments
    IDs = kwargs.get('IDs', None)
    delimiter = kwargs.get('delimiter', ',')
    # Read in the similarity scores
    sim_list = []
    with open(sim_scores) as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            sim_list.append(row)
    sim_np = np.asarray(sim_list)
    # Replace NaN scores with 0.0 if they exist.
    arr_cleaned = clean_NaN(sim_np)
    
    # Read in the IDs for ligands if given.
    if IDs is not None:
        with open(IDs, 'r') as f:
            ID_list = f.read().splitlines()
    else:
        ID_list = list(range(sim_np.shape[0]))
    return arr_cleaned, ID_list
    
    
def db_mol_IDs(db_mol, n_arr):
    '''
    Generates ID lists from lomap.DBMolecules
    for clustering or other uses.
    
        Parameters:
            db_mol: object created by lomap.DBMolecules
            n_arr: chemical similarity 2D numpy array
        
        Returns:
            ID_list: list of ID names
    '''
    mol_names = ['ID']
    for i in range(n_arr.shape[1]):
        fname = db_mol[i].getName()
        fsplit = fname.split(".", 1)[0]
        mol_names.append(fsplit)
    # Convert the name list to an np array.
    mol_arr_pre = np.array(list(mol_names))
    # Cleave the 'ID' entry
    ID_list = mol_arr_pre[1:]
    return ID_list
    

def rand_sim_scores(N):
    '''
    Random weight generator.  Generates a NxN symmetric similarity array.
    
        Parameters:
            N: the number of rows or columns
        Returns:
            b_symm: symmetric array
    '''
    b = np.random.random_integers(0,5,size=(N,N))
    # Symmetrize
    b_symm = (b + b.T)/2
    return b_symm



def clean_up(files):
    '''
    Removes leftover files.
    
        Parameters:
            files: list of files to remove.
    '''
    print("Cleaning up files from ligand posing...")
    for i in files:
        if os.path.exists(i):
            os.remove(i)
            print(f"{i} was removed.")
        else:
            print(f"The file {i} does not exist for clean up.")


class Logger(object):
    '''To use Logger, call sys.stdout = Logger().'''
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("optimal_lomap.log", "a")

    def write(self, outputs):
        self.terminal.write(outputs)
        self.log.write(outputs)

    def flush(self):
        pass
  
  
def write_json(dict, ofile):
    '''
    Write json file directly from dictionary.
    
        Parameters:
            dict: dictionary
            ofile: name of output file
    '''
    json_str = json.dumps(dict, indent=4)
    with open(ofile, 'w') as outfile:
        outfile.write(json_str)
        

def record_dicts(sub_ID, **kwargs):
    '''
    Print out the dictionaries of ligands clustered and the
    reference ligands if they were provided.
        
        Parameters:
            sub_ID: the determined clusters, IDs recorded
        
        Optional Parameters:
            sub_refs: the user entered reference ligands,
                      clustered.
               
        Output:
            Prints cluster information and saves to json
            files for future use.
    '''
    # Handle kwarg
    sub_refs = kwargs.get('sub_refs', None)
    
    print("The clusters found are:")
    # Write json file for clusters of ligands
    print(json.dumps(sub_ID, indent=4))
    write_json(sub_ID, 'cluster_IDs.json')
    # Write json file for reference ligands if indicated.
    if sub_refs is not None:
        write_json(sub_refs, 'cluster_ref_ligs.json')
        print("The passed reference ligands are in clusters:")
        print(json.dumps(sub_refs, indent=4))
    

def write_csv(data):
    '''
    Write the edge connections to a csv file.
    
        Parameters:
            data: data structure to write to csv file.
    '''
    with open('edge_data.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
