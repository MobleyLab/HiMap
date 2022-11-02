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
import lomap

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_pdf
import matplotlib.ticker as tick

from scipy import interpolate
from kneed import DataGenerator, KneeLocator

# Note: cluster_auto() is in development. Currently not exposed.

__all__ = ['dbscan', 'plt_heatmap', 'plt_dbscan', 'plt_cluster_regions',
           'plt_cluster', 'cluster_interactive', 'clusters2optimize',
           'sub_arrays']

# ****************
# FUNCTIONS
# ****************

def k_dist(data):
    '''
    Generates data for k-distance sorting to
    calculate the neighbor distance threshold for clustering.
    
        Parameters:
            data: distance array
        Returns:
            x: array the indexes each ligand
            distances: sorted nearest neighbor distances
    '''
    distances = np.sort(data, axis=0)
    distances = distances[:,1]
    x = np.linspace(0, distances.shape[0]-1, distances.shape[0])
    return x, distances


def find_shape(x, y):
    """
    Detect the direction and curvature of k-dist line.
        
        Parameters:
            x, y: x and y values
        Returns:
            direction: "increasing" or "decreasing"
            curve type: "concave" or "convex"
    """
    p = np.polyfit(x, y, deg=1)
    x1, x2 = int(len(x) * 0.2), int(len(x) * 0.8)
    q = np.mean(y[x1:x2]) - np.mean(x[x1:x2] * p[0] + p[1])
    if p[0] > 0 and q > 0:
        return 'increasing', 'concave'
    if p[0] > 0 and q <= 0:
        return 'increasing', 'convex'
    if p[0] <= 0 and q > 0:
        return 'decreasing', 'concave'
    else:
        return 'decreasing', 'convex'
 
 
def output_slopes(k, point):
    '''
    Find the slope of a curve at a x value point
        
        Parameters:
            k: in this case k will be k_raw or k_fit
            point: x point to evaluate at
            
        Returns:
            slope: slope at point on elbow/knee plots
            dslope: slope of difference curve at point on
                elblow/knee plots
    '''
    # k will be k_raw or k_fit
    x_dat, y_dat = k.x_normalized, k.y_normalized
    dx_dat, dy_dat = k.x_difference, k.y_difference
    # Find nearest to points on plots to point.
    
    def closest_vals(arr, point):
        # Find the values in array nearest to the point
        # The elif statement at 136 needs to be uncommented at 0.0 currently
        res = arr[0]
        N = len(arr)
        # Traverse the array
        for i in range(0, N, 1):
            # Find element nearest point
            if (abs(point - res) >
                abs(point - arr[i])):
                res = arr[i]
            #elif (point == arr[i]):
                #res = arr[i]
                # Find the nearest element greater
                # or less than the point.
                if res <= point:
                    try:
                        res_2 = arr[i+1]
                        index2 = i+1
                    except:
                        res_2 = arr[i-1]
                        index2 = i-1
                if res > point:
                    try:
                        res_2 = arr[i-1]
                        index2 = i-1
                    except:
                        res_2 = arr[i+1]
                        inex2 = i+1
                index = i
        # Order elements
        val1, val2 = min(res, res_2), max(res, res_2)
        i_1, i_2 = min(index, index2), max(index, index2)
        # Return elements + and - point ind indices
        return val1, val2, i_1, i_2
    
    # Collect data to find slope at point on origional curve
    x1, x2, index, index2 = closest_vals(k.x_normalized, point)
    x = [x1, x2]
    y = [k.y_normalized[index], k.y_normalized[index2]]
    eps_at_point = np.average(y)
    # Collect data to find 2nd dir at point on difference curve
    dx = [k.x_difference[index], k.x_difference[index2]]
    dy = [k.y_difference[index], k.y_difference[index2]]
    # Find the splope between the nearest points.
    slope, intercept = np.polyfit(x, y, 1)
    dslope, dintercept = np.polyfit(dx, dy, 1)
    return slope, dslope, eps_at_point

    
def find_max_curvature(x, dists, **kwargs):
    '''
    Outputs the point of maximum curvature, known as the elbow
    or knee of a curve.
        
        Parameters:
            x: x-values
            dists: distances from similarity array.
            
        Optional Parameters:
            savefigs: controls if figures are saved
            verbose: ouputs added prints to screen
            
        Returns:
            e_fit: the neighbor distance cutoff, calculated
                   or entered by user.
    '''

    # Measure concave/covex or increasing/decreasing.
    dir, cv = find_shape(x, dists)
    k_raw = KneeLocator(x, dists, S=1.0, curve=cv, direction=dir)
    k_fit = KneeLocator(x, dists, S=1.0, curve=cv, direction=dir,
                        interp_method="polynomial", online = True)

    def detect_max(k, **kwargs):
        # Find max curvature point.
        e = k.knee_y
        # If starts flat, None may be detected. Fix.
        e = 0.0 if e is None else e
        epsilon = round(e, 3)
        k.plot_knee_normalized()
        return epsilon
    try:
        e_raw = detect_max(k_raw)
    except:
        # If curve is not found, set neighbor dist to 0.5, user can override
        e_raw = 0.5
    plt.title("Raw Distance Data")
    plt.suptitle("Clustering cutoff at y-value of highest curvature",
                 size = 'medium', style = 'italic')
    plt.xlabel("Normalized ligands, sorted by distance")
    plt.ylabel("Normalized nearest neighbor distance")
    if kwargs.get('verbose') is True:
   
        plt.draw()
        plt.pause(0.1)
        input("<Hit Enter>")
        plt.savefig('cutoff_raw.pdf')
        plt.close()
    
    try:
        e_fit = detect_max(k_fit)
    except:
        # If curve is not found, set neighbor dist to 0.5, user can override
        e_fit = 0.5
    plt.title("Polynomially Fit Distance Data")
    plt.suptitle("Clustering cutoff at y-value of highest curvature",
                 size = 'medium', style = 'italic')
    plt.xlabel("Normalized ligands, sorted by distance")
    plt.ylabel("Normalized nearest neighbor distance")
    if kwargs.get('verbose') is True:
        
        plt.draw()
        plt.pause(0.1)
        input("<Hit Enter>")
        plt.savefig('cutoff_fit.pdf')
        plt.close()
        
        print(f"A suggested range for neighbor distances is between {e_raw} and {e_fit}.")
        print(f"The computed, default cutoff is {e_fit}.")
    
    
    # Testing block for analysis
    '''
    points = [0.05, 0.125, 0.175]

    for point in points:
    
        slope, dslope, eps_at_point = output_slopes(k_fit, point)
        print(f'The slope at {point} is {slope} for k_fit.')
        print(f'The slope of the diff curve at {point} is {dslope} for k_fit.')
        print(f'Epsilon at {point} = {eps_at_point}')
    
        slope2, dslope2, eps_at_point = output_slopes(k_raw, point)
        print(f'The slope at {point} is {slope2} for k_raw.')
        print(f'The slope of the diff curve at {point} is {dslope2} for k_raw.')
        print(f'Epsilon at {point} = {eps_at_point}')
    '''
    
    return e_fit


def dbscan(X, **kwargs):
    '''
    Peforms clustering using DBSCAN.
    
        Parameters:
            X: 2D distance matrix from similarity.
        Optional Parameters:
            dist_cutoff: neighbor distance cutoff
                         default = None, will calculate
            min_s: minimum sample in cluster
                   default = 1
        Returns:
            labels: array of cluster numbers by ligand
            core_samples_mask: filters clusters
            n_clusters_: the number of clusters
    '''
    # Define optional arugments
    # Minimum sample in cluster
    min_s = kwargs.get('min_s', 1)
    if min_s is None:
        min_s = 1
    else:
        min_s = min_s
    # Max distances apart to be neighbors
    dist_cutoff = kwargs.get('dist_cutoff', None)
    if dist_cutoff is None:
        # If not given, calculate it. Default.
        x, dists = k_dist(X)
        dist_cutoff = find_max_curvature(x, dists)
    else:
        dist_cutoff = dist_cutoff
        
    # Find clusters.
    db = DBSCAN(eps=dist_cutoff, min_samples=min_s, metric = 'precomputed').fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    # Find number of clusters, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    # Print cluster information for user.
    print("Estimated number of clusters: %d" % n_clusters_)
    print("Estimated number of noise points: %d" % n_noise_)
    return labels, core_samples_mask, n_clusters_


def plt_heatmap(data, ID_list, ax, fig, **kwargs):
    '''
    Plots a heatmap of the chemical distances that will be clustered.
    
        Parameters:
            data: symmetric distance array (1 - similarity)
            ID_list: list of names of IDs
            
        Optional Parameters:
            cmap: color to plot with, default is matplotlib 'CMRmap'
            tick_interval: default produces 15 or fewer ticks.
            
        Returns:
            ax: the subfigure.
    '''
    # Clean ID list if needed, func in utils.
    lomap.clean_ID_list(ID_list)
    
    # Ensure input data has 0 distance on diagonal
    np.fill_diagonal(data, 0.0)
    # Get number of ligands
    N = data.shape[0]
    # Define optional arugments
    # Plot coloring
    cmap = kwargs.get('cmap', 'inferno')

    # X and Y ticks
    # Temp note: I just deleted , 10 as the default in get kwargs
    #tick_interval = kwargs.get('tick_interval')
    #if tick_interval is None:
        # The default is to have 15 or fewer ticks.
    #if N > 15:
    tick_interval = N//15 + 1
    #else:
        #tick_interval = 1
    #else:
        #tick_interval = tick_interval

    # Heatmap plot
    im = ax.imshow(data, origin='upper', cmap='CMRmap')
    
    # Add residue ID labels to axes.
    ax.set_yticks(np.arange(N)[::tick_interval])
    ax.set_xticks(np.arange(N)[::tick_interval])
    ax.set_yticklabels(ID_list[::tick_interval])
    ax.set_xticklabels(ID_list[::tick_interval], rotation = 45, ha='right', rotation_mode='anchor')

    # Add figure labels and titles.
    plt.ylabel('Ligand IDs')
    plt.xlabel('Ligand IDs')
    plt.title('Distances in Similarity', pad=20)

    # Colorbar controls.
    im.set_clim(0, 1)
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel('Distance')
    return ax


def plt_dbscan(data, labels, core_samples_mask, n_clusters_):
    '''
    Plots clusters on 2d plot to see distance between clusters.
    
        Parameters:
            data: symmetric distance array (1 - similarity)
            labels: array of cluster numbers by ligand
            core_samples_mask: filters clusters
            n_clusters_: the number of clusters
            
        Returns:
            fig: the figure.
    '''
    fig, ax3 = plt.subplots()
    fig.set_size_inches(6, 6)
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]
        class_member_mask = labels == k
        # Plot clusters.
        xy = data[class_member_mask & core_samples_mask]
        plt.plot(
            xy[:, 0],
            xy[:, 1],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=14,
        )
        # Plot noise points.
        xy = data[class_member_mask & ~core_samples_mask]
        plt.plot(
            xy[:, 0],
            xy[:, 1],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=6,
        )
    plt.title("Estimated number of clusters: %d" % n_clusters_)
    return fig


def plt_cluster_regions(labels, ID_list, **kwargs):
    '''
    Plots the cluster regions relative to heatmap of the chemical distances.
    
        Parameters:
            labels: calculated by DBSCAN
            ID_list: list of names of IDs
            
        Optional Parameters:
            cmap: color to plot with, default is matplotlib 'CMRmap'
            tick_interval: default produces 15 or fewer ticks.
            
        Returns:
            ax: the subfigure.
    '''
    # Clean ID list if needed, func in utils.
    lomap.clean_ID_list(ID_list)
    N = len(ID_list)
    print(f'N is {N}')
    # Define optional arugments
    # Plot coloring
    cmap = kwargs.get('cmap', 'inferno')

    # Ticks parameters
    #tick_interval = kwargs.get('tick_interval')
    #if tick_interval is None:
        # The default is to have 15 or fewer ticks.
    #if N > 15:
    tick_interval = N//15 + 1
    #else:
        #tick_interval = 1
    #else:
        #tick_interval = tick_interval
    # Option to pass ax
    ax = kwargs.get('ax', None)
    if ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(4, 6)
    else:
        fig = kwargs.get('fig', None)
        ax = ax
        fig = fig
        
    # Start plotting.
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    # Make the labels a 2D array that can be plotted
    labels_arr = labels.reshape(N, 1)
    # Define discrete colors for N clusters.
    cmap = cm.get_cmap('inferno', len(unique_labels))
    # Plot. Had 0.04 for the aspect before
    psm = ax.imshow(labels_arr, cmap=cmap, rasterized=True, aspect= 'auto')
    
    # Control color bar.
    cbar = fig.colorbar(psm, ax=ax, ticks=sorted(list(unique_labels)))
    
    # What I want is to just replace -1 with noise:
    '''
    labels_l = sorted(list(unique_labels))
    for i in range(len(labels_l)):
        if labels_l[i] == -1:
             labels_l[i] = 'noise'

    if -1 in unique_labels:
        #strings = [str(x) for x in sorted(list(unique_labels))[1:]]
        cbar.ax.set_yticklabels(['noise', '0', '1', '2'])
    '''
    cbar.set_label('Cluster Number', rotation=270, labelpad = 15)
    # Add residue ID labels to axes
    ax.set_yticks(np.arange(N)[::tick_interval])
    ax.set_yticklabels(ID_list[::tick_interval])

    ax.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)

    # Add figure labels and titles
    ax.set_ylabel('Ligand IDs')
    plt.title('Cluster Regions', pad=20)
    return ax
 
 
# This need the kwargs options put in
def plt_cluster(data, labels, ID_list):
    '''
    Function to combine heatmap plot and cluster region plots into
    the same figure to compare side by side.
    
        Parameters:
            data: the distance array from similarity scores
            labels: cluster numbers calculated by DBSCAN()
            ID_list: list of ID names
             
        Returns:
            fig: the figure.
    '''
    # Clean ID list if needed, func in utils.
    lomap.clean_ID_list(ID_list)
    # This was 4.5 not 5.5
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.5), gridspec_kw={'width_ratios': [1, 4], 'height_ratios': [1]})
    fig.tight_layout()
    plt.subplots_adjust(top=0.86, bottom=0.22)
    plt_cluster_regions(labels, ID_list, ax=axs[0], fig=fig)
    plt_heatmap(data, ID_list, axs[1], fig)
    axs[0].set_title('Cluster Regions', pad = 15)
    axs[1].set_title('Heatmap of chemical distance', pad = 15)
    return fig


def clusters_w_ref(ref_ligs, sub_ID):
    '''
    Find which clusters contain reference ligands.
    Input a list of reference ligand names and check
    which cluster it is found in in dict sub_ID.
    
    Outputs: list of clusters containing a ref lig.
    
        Parameters:
            ref_ligs: list of reference ligands.
            sub_ID: dictionary, ligand IDs subdivided into clusters.
             
        Returns:
            cluster_set: list of clusters containing a ref lig.
            sub_refs: dictionary of reference ligands places by key
                      where key is the cluster number.
    '''
    hits = []
    # Make dictionary with cluster number keys to store ref ligs
    keyList = set(sub_ID.keys())
    sub_refs = {key: [] for key in keyList}
    # Find which clusters contain reference ligands.
    for i in ref_ligs:
        for k in sub_ID:
            # If ref lig name in list
            if i in sub_ID[k]:
                hits.append(k)
                # Return unique cluster keys.
                cluster_set = list(set(hits))
                # If ref in cluster, append to dict
                sub_refs[k].append(i)
    for k in keyList:
        if not sub_refs[k]:
            # Delete keys w empty lists.
            del sub_refs[k]
    return cluster_set, sub_refs
    

def cluster_auto(data, ID_list, **kwargs):
    ''' The full automated sequence, not including outputting new arrays'''
    # Clean ID list if needed, func in utils.
    lomap.clean_ID_list(ID_list)
    # Make output plots for PDF
    pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
    x, dists = k_dist(data)
    epsilon_fit = find_max_curvature(x, dists, savefigs=True)
    labels, mask, n_clusters_ = lomap.dbscan(data)
    fig1 = plt_cluster(data, labels, ID_list)
    fig2 = plt_dbscan(data, labels, mask, n_clusters_)
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.close()
    return labels
 
 
def sub_arrays(labels, n_arr, ID_list):
    '''
    Make dictionaries containing clusters to optimize. Keys are cluster
    numbers. Key contents are np similarity arrays (sub_arrs)
    or ligand names (sub_IDs).
    
        Parameters:
            labels: cluster numbers calculated by DBSCAN()
            n_arr: the 2D similarity array
            ID_list: list of ID names
             
        Returns:
            sub_arrs: dict, n_arr subdivided into dict of cluster number keys
            sub_IDs: dict, ID_list subdivided into dict of cluster number keys
    '''
    # Clean ID list if needed, func in utils.
    lomap.clean_ID_list(ID_list)
    
    # Filter out ligs measured as noise. Noise is cluster -1.
    labels_w_o_noise = [x for x in labels if x >= 0]
    unique_labels = set(labels_w_o_noise)
    
    # Generate dictionary of names for submatrices, named by cluster index (0,...N)
    sub_arrs = dict((i, "n_arr_" + str(i)) for i in range(len(unique_labels)))
    sub_IDs = dict((i, "IDs_" + str(i)) for i in range(len(unique_labels)))
    # Loop over the unique labels to generate similarity matrices of clusters
    c = 0
    for i in unique_labels:
        # Find label entries corresponding to clusters
        result = np.where(labels == i)
        # Ouput dict of clustered similarity np arrays
        sub_arrs[c] = n_arr[np.ix_(result[0][:],result[0][:])]
        # Ouput dict of clustered ID lists
        sub_IDs[c] = [ID_list[index] for index in result[0][:]]
        c = c + 1
    return sub_arrs, sub_IDs
    
 
def cluster_interactive(sim_data, ID_list):
    '''
    Function for if user wants to inspect distance data first and
    self assign the clustering neighbor distance cutoff. Also useful if
    testing different potential cutoff values. This is the current default.
    
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
    '''
    def get_numeric(splitme):
        '''
        Will read str and return integer values.
        '''
        delimiters = [" ", ",", "[", "]", ")", "("]
        cleaned = splitme
        for i in delimiters:
            cleaned = cleaned.replace(i, " ")
        digs = [int(s) for s in cleaned.split() if s.isdigit()]
        return digs
    
    header = "                    Cluster Selection table                      "
    line = "-----------------------------------------------------------------"
    
    # Clean ID list if needed, func in utils.
    lomap.clean_ID_list(ID_list)
    
    # Generate distance data
    data = 1 - sim_data
    
    # Output distance info to user.
    x, dists = k_dist(data)
    auto_cutoff = find_max_curvature(x, dists, savefigs=True, verbose=True)
    # Temp putting this here to save figs.
    #auto_cutoff = find_max_curvature(x, dists, savefigs=True, verbose=False)
    #plt.close()
    
    
    # Get user input on distance cutoff, epsilon in dbscan, typed
    input_cutoff = input("Enter a neighbor distance cutoff for clustering:")
    if input_cutoff == "":
        # If enter, set to autodetected cutoff
        user_cutoff = auto_cutoff
    else:
        user_cutoff = float(input_cutoff)
    # Do Clustering.
    # I changed the min_s from 4 to 2
    labels, mask, n_clusters_ = dbscan(data, dist_cutoff=user_cutoff, min_s = 2)
    #ax = plt_cluster_regions(labels, ID_list)
    ax = plt_cluster(data, labels, ID_list)
    # Output table for user cluster selection
    print("Which clusters do you want to optimize?")
    print(line)
    print(header)
    print(line)
    d = {'a': ["all", 'Optimize all clusters.'],
         'w': ["w_ref_lig", 'Optimize clusters w/ reference ligands.'],
         'ints': ["list of ints", 'Select cluster numbers from figure.']
        }
    print ("{:<8} {:<15} {:<10}".format('Input','Variable','Meaning'))
    for k, v in d.items():
        vari, meani = v
        print ("{:<8} {:<15} {:<10}".format(k, vari, meani))
    print(line)

    plt.tight_layout()
    plt.draw()
    plt.pause(0.1)
    input("<Hit Enter>")
    plt.close()
    user_clusters = input("Input a, w, or cluster numbers(ints):")
    
    # Figure saving
    pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
    fig1 = plt_cluster(data, labels, ID_list)
    fig2 = plt_dbscan(data, labels, mask, n_clusters_)
    pdf.savefig(fig1, bbox_inches= "tight", pad_inches=0.5)
    pdf.savefig(fig2)
    pdf.close()
    
    # Take use input and translate to kwarg option.
    if user_clusters == "" or user_clusters == 'a' or user_clusters == 'all':
        # If only enter, set to default of 'all'
        selected_clusters = 'all'
    elif user_clusters == 'w' or user_clusters == 'w_ref_lig':
        selected_clusters = 'w_ref_lig'
    elif get_numeric(user_clusters):
        selected_clusters = get_numeric(user_clusters)
        # Test if entered ints are valid.
        if set(selected_clusters).issubset(set(labels)) is False:
            non_noise = [item for item in set(labels) if item >= 0]
            raise ValueError(f"Invalid cluster numbers. Valid numbers are {non_noise}")
    else:
        raise ValueError("Enter valid input from Cluster Selection table, example: 1, 2")
    
    # Generate sub-arrays of clusters. Stored in dictionaries.
    sub_arr, sub_ID = lomap.sub_arrays(labels, sim_data, ID_list)

    return sub_arr, sub_ID, selected_clusters
 
 
def clusters2optimize(sub_arr, sub_ID, **kwargs):
    '''
    Sends clusters for optimization. Takes lomap.Optimize() kwargs.
    
        Parameters:
            sub_arr = dictionary of sub arrays of similarity for each cluster.
            output from lomap.sub_arrays()
            sub_ID = dictionary of sub arrays of ligand ID names over each cluster.
            output from lomap.sub_arrays()
            
        Optional Parameters, passed for optimization:
            clusters2optim = vector of cluster numbers to optimize.
                (ex: [1, 4]) default: 'all'
                If 'all' is chosen, all calculated clusters will be iteratively
                optimized.
                
            db_mol      = output of lomap.DBMolecules()
            
            optim_types = [ optim1 , optim2 ](str), optimization type such as
                        'A' and 'D'. Currently two types are required. Options are:
                        'A', 'D', 'P', 'mA', 'mP', 'negA', 'negD', 'random'
                        
            ref_ligs (list of str) = user input reference ligand. If not manually
                selected, will be calculated based on max similarity. Multiple can
                be accepted.
                
            ID_list (str) = the user can define ligand names. The default
                is not not enter a list but instead output.
                
            num_edges = the number of edges requested for optimization. If 'n' is the
                        number of ligands in each cluster, the options are:
                        '1n', '2n', 'nlnn', 'min', 'max', integers.
                        The edge number requested must be in the range of [min, max].
                        If less than min, will be set to min. If greater than max,
                        will be set to max.
        
        Outputs:
            Optimal graph outputs.
                
    Example usages:
    lomap.clusters2optimize(sub_arr, sub_ID, clusters2optim = 'all', optim_types = ['A', 'D'])
    lomap.clusters2optimize(sub_arr, sub_ID)
    '''
    clusters2optim = kwargs.get('clusters2optim', 'all')
    if clusters2optim is 'all':
        sub_arr_passed = sub_arr
        sub_ID_passed = sub_ID
        if 'ref_ligs' in kwargs:
            ref_ligs = kwargs.get('ref_ligs')
            # Find which ref ligs are in which cluster if passed.
            clusts_w_ref, sub_refs = clusters_w_ref(ref_ligs, sub_ID)
            
    elif clusters2optim is 'w_ref_lig':
        if 'ref_ligs' not in kwargs:
            raise ValueError("if clusters='w_ref_lif', you must input ref_ligs kwarg.")
        else:
            ref_ligs = kwargs.get('ref_ligs')
            print(ref_ligs)
            # Find which clusters contain reference ligands and ref_ligs/cluster.
            clusts_w_ref, sub_refs = clusters_w_ref(ref_ligs, sub_ID)
            clusters_input = tuple(clusts_w_ref)
            # Retain only those clusters with reference ligands.
            sub_arr_passed = {key: sub_arr[key] for key in clusters_input}
            sub_ID_passed = {key: sub_ID[key] for key in clusters_input}
    
    # If user entered list of cluster integers to optimize.
    elif isinstance(clusters2optim, list):
        # Test if passed cluster keys are allowed.
        clusters_input = tuple(clusters2optim)
        for i in clusters_input:
            if i not in sub_arr.keys():
                raise ValueError(f'Accepted cluster ints are {sub_arr.keys()}')
                
        # Extract user selected clusters from dict.
        sub_arr_passed = {key: sub_arr[key] for key in clusters_input}
        sub_ID_passed = {key: sub_ID[key] for key in clusters_input}
        # If user passed in reference ligands, extract.
        if 'ref_ligs' in kwargs:
            ref_ligs = kwargs.get('ref_ligs')
            # Find which ref ligs are in which cluster.
            clusts_w_ref, sub_refs = clusters_w_ref(ref_ligs, sub_ID)
    else:
        raise ValueError("clusters2optim = 'all', 'w_ref_lif', or [ints].")
        
    # If no ref_ligs given, send each passed cluster for optimization.
    if 'ref_ligs' not in kwargs:
        # Output jsons of cluster info
        lomap.record_dicts(sub_ID)
        for i,j in zip(sub_arr_passed, sub_ID_passed):
            n_ar = sub_arr_passed[i]
            sub_ID_list = sub_ID_passed[j]
            lomap.Optimize(n_ar, ID_list = sub_ID_list, **kwargs)
    else:
        # There are ref_ligs in some or all clusters.
        # If clusters2optim = 'all' or [ints], passed clusters may have ref ligs.
        # If clusters2optim = 'w_ref_lig', each passed cluster will have a ref lig.
        
        # Output jsons of cluster info
        lomap.record_dicts(sub_ID, sub_refs = sub_refs)
        
        # This should always be true
        assert set(clusts_w_ref).issubset(sub_refs.keys())
        # Optimize clusters with and without reference ligands found in cluster.
        for k in sub_ID_passed.keys():
            n_ar = sub_arr_passed[k]
            sub_ID_list = sub_ID_passed[k]
            if k in set(clusts_w_ref):
                # If the cluster has a reference ligand, optimize with ref_lig.
                lomap.Optimize(n_ar, ID_list = sub_ID_list, ref_lig = sub_refs[k], **kwargs)
            else:
                # If the cluster has no reference ligand, optimize.
                lomap.Optimize(n_ar, ID_list = sub_ID_list, **kwargs)
