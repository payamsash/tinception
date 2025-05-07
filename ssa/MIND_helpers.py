from scipy.spatial import cKDTree as KDTree
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict
from itertools import combinations

def is_outlier(points, thresh=7): #taken from https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data

    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def get_KDTree(x): #Inspired by https://gist.github.com/atabakd/ed0f7581f8510c8587bc2f41a094b518

    # Check the dimensions are consistent
    x = np.atleast_2d(x)
    
    # Build a KD tree representation of the samples
    xtree = KDTree(x)
    
    return xtree

def get_KL(x, y, xtree, ytree): #Inspired by https://gist.github.com/atabakd/ed0f7581f8510c8587bc2f41a094b518

    
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)

    n,d = x.shape
    m,dy = y.shape
    
    #Check dimensions
    assert(d == dy)

    # Get the first two nearest neighbours for x, since the closest one is the
    # sample itself.
    r = xtree.query(x, k=2, eps=.01, p=2)[0][:,1]
    s = ytree.query(x, k=1, eps=.01, p=2)[0]
    
    rs_ratio = r/s

    #Remove points with zero, nan, or infinity. This happens when two regions have a vertex with the exact same value – an occurence that basically onnly happens for the single feature MSNs
    #and has to do with FreeSurfer occasionally outputting the exact same value for different vertices.
    rs_ratio = rs_ratio[np.isfinite(rs_ratio)]
    rs_ratio = rs_ratio[rs_ratio!=0.0]
    
    # There is a mistake in the paper. In Eq. 14, the right side misses a negative sign
    # on the first term of the right hand side.

    kl = -np.log(rs_ratio).sum() * d / n + np.log(m / (n - 1.))
    kl = np.maximum(kl, 0)
    
    return kl


def calculate_mind_network(data_df, feature_cols, region_list, resample=False, n_samples = 4000):

    '''
    MIND = pd.DataFrame(np.zeros((len(region_list), len(region_list))), \
                        index = region_list, columns = region_list)

    #Get only desired regions
    data_df = data_df.loc[data_df['Label'].isin(region_list)]
    
    #Resample dataset if resample has been set to True and if it is UNIVARIATE ONLY. This should only be done if you are using a single feature which contains repeated values.
    if (len(feature_cols) == 1) and resample==True:   
        n_samples = n_samples
        resampled_dataset = pd.DataFrame(np.zeros((n_samples, len(region_list))), columns = region_list)

        for name, data in data_df.groupby('Label'):
            resampled_dataset[name] = stats.gaussian_kde(data[feature_cols[0]]).resample(n_samples)[0]

        resampled_dataset = resampled_dataset.melt(var_name = 'Label', value_name = feature_cols[0])
        data_df = resampled_dataset

    if (len(feature_cols) > 1) and resample==True:   
        raise Exception("Resampling the data is only supported if you are using a single feature -- this is because higher order density estimation can be unreliable and very computationally expensive.")

    #Check that there aren't many repeated values
    percent_unique_vals = len(data_df[feature_cols].drop_duplicates())/len(data_df[feature_cols])
    
    if percent_unique_vals < 0.8:
        raise Exception("There are many repeated values in the data, which compromises the validity of MIND calculation. Please minimize the number of repeated values in the data and try again. If you are using only one feature, try rerunning with resample=True.")

    grouped_data = data_df.groupby('Label')

    KDtrees = defaultdict(object)

    for i, (name_x, dat_x) in enumerate(grouped_data):
        tree = get_KDTree(dat_x[feature_cols])
        KDtrees[name_x] = tree
    
    used_pairs = []
    
    for i, (name_x, dat_x) in enumerate(grouped_data):
        
        for name_y, dat_y in grouped_data:
            if name_x == name_y:
                continue

            if set([name_x,name_y]) in used_pairs:
                continue

            dat_x = dat_x[feature_cols]
            dat_y = dat_y[feature_cols]

            KLa = get_KL(dat_x, dat_y, KDtrees[name_x], KDtrees[name_y])
            KLb = get_KL(dat_y, dat_x, KDtrees[name_y], KDtrees[name_x])

            kl = KLa + KLb

            MIND.at[name_x,name_y] = 1/(1+kl)
            MIND.at[name_y,name_x] = 1/(1+kl)

            used_pairs.append(set([name_x,name_y]))

    MIND = MIND[region_list].T[region_list].T
    
    return MIND
    '''

    MIND = pd.DataFrame(np.zeros((len(region_list), len(region_list))), 
                    index=region_list, columns=region_list)

    # Filter once
    data_df = data_df[data_df['Label'].isin(region_list)]

    # Group data and build KDTree for each region
    grouped_data = {label: df[feature_cols].values for label, df in data_df.groupby('Label')}
    KDtrees = {label: get_KDTree(features) for label, features in grouped_data.items()}

    # Iterate only over unique unordered region pairs
    for name_x, name_y in combinations(region_list, 2):
        dat_x = grouped_data[name_x]
        dat_y = grouped_data[name_y]

        KLa = get_KL(dat_x, dat_y, KDtrees[name_x], KDtrees[name_y])
        KLb = get_KL(dat_y, dat_x, KDtrees[name_y], KDtrees[name_x])

        kl = KLa + KLb
        sim = 1 / (1 + kl)

        MIND.at[name_x, name_y] = sim
        MIND.at[name_y, name_x] = sim

    # Optional: ensure final MIND matrix is properly ordered
    MIND = MIND.loc[region_list, region_list]

    return MIND