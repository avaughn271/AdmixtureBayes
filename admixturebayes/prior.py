from scipy.stats import geom, nbinom
from Rtree_operations import (get_all_branch_lengths, get_all_admixture_proportions, get_number_of_leaves)
from math import log
from uniform_topological_prior import uniform_topological_prior_function
from Rtree_to_covariance_matrix import get_admixtured_populations

def logpdf(x, fro=0.0, to=1.0):
    log_max_height=log((to-fro)*2)
    log_frac=log( (to-(x-fro))/(to-fro) )
    return log_max_height+log_frac

def calculate_branch_prior(branches, n):
    rate=float(2*n-2)/len(branches)
    n2k=len(branches)
    d=float(n2k)/float(n)
    
    return -sum(branches)*rate+log(rate)*len(branches)

def prior(x, p=0.5, pks={}, r=0):
    tree, add=x
    no_leaves=get_number_of_leaves(tree)
    admixtures=get_all_admixture_proportions(tree)
    if not all(prop>=0 and prop <=1 for prop in admixtures):
        return -float('inf')
    branches=get_all_branch_lengths(tree)
    if not all(branch>=0 for branch in branches):
        return -float('inf')
    branch_prior=calculate_branch_prior(branches, no_leaves)
    no_admix_prior=no_admixes(p, len(admixtures), r=r)
    admix_prop_prior=0
    top_prior=uniform_topological_prior_function(tree)
    logsum=branch_prior+no_admix_prior+admix_prop_prior+top_prior-add
    pks['branch_prior']= branch_prior
    pks['no_admix_prior']=no_admix_prior
    pks['admix_prop_prior']=admix_prop_prior
    pks['top_prior']= top_prior
    return logsum

def linear_admixture_proportions(admixtures):
    return sum((logpdf(admixture) for admixture in admixtures))

def no_admixes(p, admixes, hard_cutoff=20, r=0):
    if admixes>hard_cutoff:
        return -float('inf')
    if r>1:
        if hard_cutoff is None:
            return nbinom.logpmf(admixes,n=r, p=1.0-p)
        else:
            return nbinom.logpmf(admixes, n=r, p=1.0 - p) - nbinom.logcdf(hard_cutoff ,n=r, p= 1.0 - p)
    else:
        if hard_cutoff is None:
            return geom.logpmf(admixes+1, 1.0-p)

        return geom.logpmf(admixes+1, 1.0-p)-geom.logcdf(hard_cutoff+1, 1.0-p)

def matchmake(single_coalescences, coalescences_on_hold):
    happy_couples=[]
    continuing_singles=[]
    for key,branch in coalescences_on_hold:
        if (key,branch) in single_coalescences:
            happy_couples.append(((key,branch),single_coalescences[(key,branch)]))
            del single_coalescences[(key,branch)]
        else:
            continuing_singles.append((key,branch))
    return single_coalescences, happy_couples, continuing_singles