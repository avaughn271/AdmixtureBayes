from Rtree_operations import get_trivial_nodes
from reduce_covariance import reduce_covariance, thin_covariance

from numpy import loadtxt, array, insert, amin, delete, ix_, ones,nan, dtype
from copy import deepcopy
import subprocess
import os

def read_freqs(new_filename, locus_filter):
    with open(new_filename, 'r') as f:
        names=f.readline().split()
        allele_counts=[]
        pop_sizes=[]
        minors=[]
        total_sum=0
        for n,r in enumerate(f.readlines()):
            minor_majors=r.split()
            minor_list=[]
            freqs=[]
            pop_sizes_SNP=[]
            for minor_major in minor_majors:
                minor, major= list(map(float,minor_major.split(',')))
                total_sum+=major+minor
                if major+minor==0:
                    freqs.append(nan)
                else:
                    freqs.append(float(minor)/float(major+minor))
                minor_list.append(minor)
                pop_sizes_SNP.append(major+minor)
            if locus_filter(freqs,pop_sizes, names):
                minors.append(minor_list)
                pop_sizes.append(pop_sizes_SNP)
                allele_counts.append(freqs)
    return allele_counts, names, pop_sizes, minors, total_sum


def get_xs_and_ns_from_freqs(ps, npop, locus_filter):
    mat=[]
    names=[]
    for k,p in list(ps.items()):
        mat.append(p)
        names.append(k)
    mat=array(mat)
    ns=npop*ones(mat.shape)
    mat, ns= locus_filter.apply_filter(mat, ns, names)
    xs=mat*npop
    return xs,ns,names

def get_xs_and_ns_from_treemix_file(snp_file, locus_filter):
    if snp_file.endswith('.gz'):
        new_filename=unzip(snp_file)
    else:
        new_filename=snp_file
    allele_freqs, names, ns, minors, total_sum= read_freqs(new_filename, locus_filter)
    xs=array(minors, dtype=dtype(float)).T
    ns=array(ns, dtype=dtype(float)).T
    return xs,ns,names

def order_covariance(xnn_tuple, outgroup=''):
    if not outgroup:
        return xnn_tuple
    xs,ns,names=xnn_tuple
    assert outgroup in names, 'The outgroup was not found in the data. Did you spell it correctly?'
    n_outgroup=next((n for n, e in enumerate(names) if e==outgroup))
    xs_o=xs[n_outgroup,:]
    ns_o=ns[n_outgroup,:]
    names_o=names[n_outgroup]
    xs=delete(xs, n_outgroup,0)
    ns=delete(ns, n_outgroup,0)
    names.remove(names_o)
    xs=insert(xs, 0, xs_o, axis=0)
    ns=insert(ns, 0, ns_o, axis=0)
    names=[names_o]+names
    return xs,ns,names

def _get_permutation(actual, target):
    val_to_index={val:key for key, val in enumerate(actual)}
    indices_to_get_target=[val_to_index[val] for val in target]
    return indices_to_get_target

def reorder_covariance(cov, names, full_nodes):
    indices=_get_permutation(names, full_nodes)
    return cov[ix_(indices,indices)]

def reorder_reduced_covariance(cov, names, full_nodes, outgroup=''):
    n1=len(names)
    n2=len(full_nodes)
    m=cov.shape[0]
    assert n1==n2==(m+1), 'Unexpected input'
    names2=deepcopy(names)
    full_nodes2=deepcopy(full_nodes)
    names2.remove(outgroup)
    full_nodes2.remove(outgroup)
    indices=_get_permutation(names2, full_nodes2)
    return cov[ix_(indices, indices)]

def normalise(m):
    return m-max(0,amin(m))

def file_to_emp_cov(filename, reduce_column=None, nodes=None, sort_nodes_alphabetically=False, vc=None, return_only_covariance=True, subnodes=[]):
    dat=[]
    multiplier=None
    with open(filename, 'r') as f:
        actual_nodes=f.readline().rstrip().split(" ")
        for i,l in enumerate(f.readlines()):
            if i>=len(actual_nodes):
                if len(l)>4:
                    multiplier=float(l.split('=')[1])
                break
            n=list(map(float, l.split()[1:]))
            if len(n)>1:
                dat.append(n)

    m=array(dat)
    if vc:
        varc=loadtxt(vc)
    mapping={val:key for key, val in enumerate(actual_nodes)}
    #removedprin 'nodes', nodes
    if nodes is None and sort_nodes_alphabetically:
        nodes=sorted(actual_nodes)
    if nodes is not None:
        new_order=[mapping[node] for node in nodes]
        #removedprin 'm.shape', m.shape
        m=m[:, new_order][new_order]
        if vc:
            varc=varc[:,new_order][new_order]
    if reduce_column is not None:
        m=reduce_covariance(m, reduce_column)
        #m=normalise(m)
    if subnodes:
        m=thin_covariance(m, nodes, subnodes)
        if vc:
            varc=thin_covariance(vc, nodes,subnodes)
    res=[m]
    if multiplier is not None:
        res.append(multiplier)
        if vc:
            res.append(varc)
        else:
            res.append(None)

    if len(res)==1 or return_only_covariance:
        return m
    else:
        return res

def emp_cov_to_file(m, filename='emp_covimport', nodes=None):
    if nodes is None:
        n=m.shape[0]
        nodes=get_trivial_nodes(n)
    with open(filename, 'w') as f:
        f.write(' '.join(nodes)+'\n')
        for i, node in enumerate(nodes):
            f.write(node+ ' '+ ' '.join(map(str, m[i]))+'\n')
    #removedprin 'wrote matrix to file', filename

def unzip(filename, overwrite=False, new_filename=None):
    original_filename=filename[:]
    assert filename.endswith('.gz'), 'file with non-zipped ending was passed to the unzip function'
    if new_filename is None:
        new_filename=filename[:-3]
    if (not overwrite) and os.path.exists(new_filename):
        return new_filename
    command=['gunzip','-c',filename]
    #removedprin command
    with open(new_filename, 'w') as f:
        subprocess.call(command, stdout=f)
    return new_filename

def gzip(filename, overwrite=False, new_filename=None):
    original_filename = filename[:]
    assert not filename.endswith('.gz'), 'file with zipped ending was passed to the zip function'
    if new_filename is None:
        new_filename=filename+'.gz'
    if (not overwrite) and os.path.exists(new_filename):
        return new_filename
    command=['gzip','-c',filename]
    #removedprin command
    with open(new_filename, 'w') as f:
        subprocess.call(command, stdout=f)
    return new_filename
