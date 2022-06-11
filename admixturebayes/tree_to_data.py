from Rtree_to_covariance_matrix import make_covariance
from Rtree_operations import (find_rooted_nodes, get_number_of_leaves, get_real_parents, get_no_leaves,
                               node_is_admixture, node_is_leaf_node, node_is_coalescence,
                              get_real_children_root, get_trivial_nodes, scale_tree_copy, get_leaf_keys,
                               get_max_timing, get_all_branch_lengths)
from reduce_covariance import reduce_covariance, thin_covariance

from numpy import loadtxt, array, mean, vstack, sum, insert, hstack, vsplit, amin, delete, ix_, ones,nan, dtype
from copy import deepcopy
import subprocess
from scipy.stats import wishart
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

def tree_to_data_perfect_model(tree, df):
    m=make_covariance(tree)
    r=m.shape[0]
    m=wishart.rvs(df=r*df-1, scale=m/(r*df))
    return m

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
            #removedprin l
            n=list(map(float, l.split()[1:]))
            if len(n)>1:
                dat.append(n)

    m=array(dat)
    if vc:
        varc=loadtxt(vc)
    #removedprin m
    mapping={val:key for key, val in enumerate(actual_nodes)}
    #removedprin 'mapping', mapping
    #removedprin 'nodes', nodes
    if nodes is None and sort_nodes_alphabetically:
        nodes=sorted(actual_nodes)
    if nodes is not None:
        new_order=[mapping[node] for node in nodes]
        #removedprin 'new_order', new_order
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

def scaled_tupled_branches(tree, d):
    '''
    g
    '''
    for key in list(tree.keys()):
        tree[key][3]=(tree[key][3][0], tree[key][3][1]*d)
        if tree[key][4] is not None:
            tree[key][4]=(tree[key][4][0], tree[key][4][1]*d)
    return tree


def extend_branch_lengths(tree, times):
    for key, node in list(tree.items()):
        for n,parent_key in enumerate(get_real_parents(node)):
            pseudo_time=times[parent_key]-times[key]
            tree[key][3+n]=(tree[key][3+n], pseudo_time)
    return tree

def construct_ej_es_string(tree, times, leaf_keys, final_pop_size=1.0):
    s_times=sorted([(v,k) for k,v in list(times.items())])
    dic_of_lineages={(key,0):(n+1) for n,key in enumerate(leaf_keys)}
    #removedprin dic_of_lineages
    population_count=len(dic_of_lineages)
    res_string=''
    for time,key in s_times:
        if key=='r':
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(final_pop_size)+' '
            break
        node=tree[key]
        if node_is_coalescence(node):
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
        if node_is_admixture(node):
            population_count+=1
            i=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))[0]
            res_string+='-es '+str(time)+' '+str(i)+' '+str(1.0-node[2])+' '
            dic_of_lineages[(key,0)]=i
            dic_of_lineages[(key,1)]=population_count
    return res_string


def construct_ej_en_es_string(tree, times, leaf_keys, final_pop_size=1.0):
    s_times=sorted([(v,k) for k,v in list(times.items())])
    dic_of_lineages={(key,0):(n+1) for n,key in enumerate(leaf_keys)}
    #removedprin dic_of_lineages
    population_count=len(dic_of_lineages)
    res_string=''
    add_on=1e-6
    for time,key in s_times:
        if key=='r':
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            pop_size=final_pop_size
            res_string+='-en '+str(time+add_on)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
            break
        node=tree[key]
        if node_is_leaf_node(node):
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_coalescence(node):
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time+add_on)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_admixture(node):
            population_count+=1
            i=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))[0]
            res_string+='-es '+str(time)+' '+str(i)+' '+str(node[2])+' '# WRONG EARLIER VERSION: str(1.0-node[2])+' '
            dic_of_lineages[(key,0)]=i
            dic_of_lineages[(key,1)]=population_count
            pop_size1=calculate_pop_size(node[3])
            pop_size2=calculate_pop_size(node[4])
            res_string+='-en '+str(time+add_on)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size1)+' '
            res_string+='-en '+str(time+add_on)+' '+str(dic_of_lineages[(key,1)])+' '+str(pop_size2)+' '
    return res_string


def get_affected_populations(dic_of_lineages, children_branches):
    return [dic_of_lineages[children_branch] for children_branch in children_branches]

def calculate_pop_size(tup):
    drift, actual=tup
    return actual/drift*2

def call_ms_string(ms_string, sequence_file):
    with open(sequence_file, 'w') as f:
        print('running', ms_string)
        p = subprocess.Popen(ms_string, stdout=subprocess.PIPE, shell=True)
        line_number = 0
        for line in p.stdout.readlines():
            line_number += 1
            if line_number >= 5 and line and (line[0]=='0' or line[0]=='1'):
                #removedprin len(line)
                #removedprin line_number, line[:4]
                f.write(line)
            #else:
                #removedprin line_number,':', line.rstrip()
    return 0

def unzip(filename, overwrite=False, new_filename=None):
    original_filename=filename[:]
    assert filename.endswith('.gz'), 'file with non-zipped ending was passed to the unzip function'
    if new_filename is None:
        new_filename=filename[:-3]
    if (not overwrite) and os.path.exists(new_filename):
        #warnings.warn('Not unzipping '+original_filename + ' because '+ new_filename+ ' already exists')
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
        #warnings.warn('Not zipping ' + original_filename + ' because ' + new_filename + ' already exists')
        return new_filename
    command=['gzip','-c',filename]
    #removedprin command
    with open(new_filename, 'w') as f:
        subprocess.call(command, stdout=f)
    return new_filename
