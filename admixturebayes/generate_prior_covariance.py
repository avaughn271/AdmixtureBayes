from generate_prior_trees import generate_phylogeny
from Rtree_to_covariance_matrix import make_covariance
from scipy.stats import uniform, beta
from Rtree_operations import scale_tree

def generate_covariance(size, scale_metod='beta', return_tree=False):
    tree= generate_phylogeny(size)
    cov=make_covariance(tree)
    s=calc_s(scale_metod)
    if return_tree:
        return cov*s, scale_tree(tree, s)
    return cov*s

def calc_s(scale_method):
    if scale_method=='uniform':
        return uniform.rvs()
    elif scale_method=='beta':
        return beta.rvs(a=1, b=5)
    elif scale_method=='None':
        return 1
    else:
        assert False, 'unexpected scale_method: '+str(scale_method)
    