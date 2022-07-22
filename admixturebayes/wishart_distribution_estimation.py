from numpy.random import choice
from tree_to_data import unzip, gzip
import warnings

from construct_covariance_choices import empirical_covariance_wrapper_directly
from pathos.multiprocessing import Pool
import os

import numpy as np

#maximizes the function, function in the interval [lower_limit,\infty).
def I_cant_believe_I_have_to_write_this_function_myself(function, lower_limit):
    old_x=lower_limit
    new_x=lower_limit*2
    old_y=function(old_x)
    lower_y=old_y
    new_y=function(new_x)
    sgn=1
    max_step_size_increases=20
    step_size_increases=0
    step_size=lower_limit
    for _ in range(20):
        c=0
        while new_y<old_y and c<100:
            old_x=new_x
            new_x+=sgn*step_size
            old_y=new_y
            if new_x<lower_limit:
                new_x=lower_limit
                new_y=lower_y
                break
            new_y=function(new_x)
            c+=1
            if c>10 and max_step_size_increases>step_size_increases:
                c=0
                step_size*=2
                step_size_increases+=1
        sgn=-sgn
        step_size*=0.5
        old_x=new_x
        new_x=old_x+sgn*step_size
        old_y=new_y
        new_y=function(new_x)
    return new_x

def variance_mean_based(sample_of_matrices, divisor=None, verbose_level='normal'):
    mean_wishart=  np.mean(sample_of_matrices, axis=0)
    var_wishart= np.var(sample_of_matrices, axis=0)
    r=mean_wishart.shape[0]
    var_rom_mean_wishart=np.square(mean_wishart)+np.outer(np.diag(mean_wishart),np.diag(mean_wishart))
    def penalty_function(df_l):
        df=df_l
        val=np.linalg.norm(var_wishart-var_rom_mean_wishart/df)**2
        return np.log(val)
    
    rval=I_cant_believe_I_have_to_write_this_function_myself(penalty_function, r)
    return rval

    
def get_partitions(lines, blocksize):
    list_of_lists=[]
    for i in range(0,len(lines)-blocksize, blocksize):
        list_of_lists.append(lines[i:(i+blocksize)])
    if len(lines)%blocksize == 0:
        list_of_lists.append(lines[-blocksize:])
    return list_of_lists

def bootstrap_indices(k):
    bootstrap_inds=choice(k, size=k, replace = True)
    return bootstrap_inds

def combine_covs(tuple_covs, indices):
    cov_sum, scale_sum=0.0,0.0
    for i in indices:
        cov_sum+=tuple_covs[i][0]
        scale_sum+=tuple_covs[i][1]
    return cov_sum/scale_sum

def bootsrap_combine_covs(covs, cores, bootstrap_samples):
    p=Pool(cores)
    def t(empty):
        indices=bootstrap_indices(len(covs))
        return combine_covs(covs, indices)
    result_covs=list(map(t, list(range(bootstrap_samples))))
    return result_covs

def remove_files(filenames):
    for fil in filenames:
        os.remove(fil)
        os.remove(fil[:-3])

def make_covariances(filenames, cores, **kwargs):
    covs=[]
    p=Pool(cores)
    def t(filename):
        return empirical_covariance_wrapper_directly(filename, **kwargs)
    covs=list(map(t, filenames))
    try:
        remove_files(filenames)
    except OSError as e:
        warnings.warn('Erasing the files did not succeed',UserWarning)
    return covs

def make_single_files(filename,blocksize, no_blocks, prefix='', verbose_level='normal'):
    assert (blocksize is not None) or (no_blocks is not None), 'Has to specify either block size or number of blocks'
    filenames=[]
    os.mkdir(os.getcwd() + "/temp_adbayes")
    prefix = os.getcwd() + "/temp_adbayes/"
    if filename.endswith('.gz'):
        filename=unzip(filename)
    filename_reduced=prefix+filename.split(os.sep)[-1]+'boot.'
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
    if no_blocks is not None:
        blocksize=n/no_blocks
    line_sets=get_partitions(lines, blocksize)
    if verbose_level!='silent':
        print('total number of SNPs: '+ str(n))
    for i,lins in enumerate(line_sets):
        new_filename= filename_reduced+str(i)
        with open(new_filename, 'w') as g:
            g.write(first_line)
            g.writelines(lins)
        gzipped_filename=gzip(new_filename, overwrite=True)
        filenames.append(gzipped_filename)
    return filenames, first_line.split()

def estimate_degrees_of_freedom_scaled_fast(filename,
                                            bootstrap_blocksize=1000,
                                            no_blocks=None,
                                            no_bootstrap_samples=10,
                                            cores=1,
                                            prefix='',
                                            verbose_level='normal',
                                            **kwargs):
    single_files, nodes=make_single_files(filename, blocksize=bootstrap_blocksize, no_blocks=no_blocks, prefix=prefix, verbose_level=verbose_level)
    assert len(single_files)>1, 'There are ' +str(len(single_files)) + ' bootstrapped SNP blocks and that is not enough. Either add more data or lower the --bootstrap_blocksize'
    if len(single_files)<39:
        warnings.warn('There are only '+str(len(single_files))+' bootstrap blocks. Consider lowering the --bootstrap_blocksize or add more data.', UserWarning)
    single_covs=make_covariances(single_files, cores=cores, return_also_mscale=True, **kwargs)
    covs=bootsrap_combine_covs(single_covs, cores=cores, bootstrap_samples=no_bootstrap_samples)
    res=variance_mean_based(covs, verbose_level=verbose_level)
    return res, covs
