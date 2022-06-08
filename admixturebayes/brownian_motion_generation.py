from Rtree_operations import (find_rooted_nodes, node_is_admixture, node_is_leaf_node, 
                              get_branch_length, get_admixture_proportion_from_key, get_children,
                              mother_or_father, other_branch)
from scipy.stats import uniform, norm, binom, multivariate_normal
from copy import deepcopy
from numpy import clip, cov, array, mean, sqrt, zeros

def simulate_xs_and_ns(n,N, Sigma, ns, normal_xval=False):
    p0s=uniform.rvs(size=N)
    pij=zeros((n+1,N))
    for s,p0 in enumerate(p0s):
        pij[1:,s]=multivariate_normal.rvs(mean=[p0]*n, cov=Sigma*p0*(1-p0))
        pij[0,s]=p0
    pij2=clip(pij,0,1)
    ##removedprin ns
   # print pij
    if normal_xval:
        xs=norm.rvs(loc= ns*pij2, scale=np.sqrt(pij2*(1-pij2)*ns) ) 
    else:
        xs=binom.rvs(ns.astype(int), p=pij2)
    #xs=np.clip(xs, 0,ns)
    return xs, p0s, pij

def add_noise(p,branch_length):
    #return add_noise2(p, branch_length)
    p=array(p)
    sd=sqrt(branch_length*p*(1-p))
    p=[a+norm.rvs(scale=sd[n]) if 0<a<1 else a for n,a in enumerate(p) ]
    p=clip(p,0,1)
    return p

def add_noise3(p,branch_length):
    #return add_noise2(p, branch_length)
    sd=sqrt(branch_length)
    p=[a+norm.rvs(scale=sd) if 0<a<1 else a for n,a in enumerate(p) ]
    p=clip(p,0,1)
    return p

def remove_non_snps(p, outgroup=''):
    names=[]
    vals=[]
    for k,v in list(p.items()):
        names.append(k)
        vals.append(v)
    print(names)
    if outgroup:
        n_outgroup=next((n for n, e in enumerate(names) if e==outgroup))
        print(n_outgroup)
        ad=vals[n_outgroup]                                  
    else:
        ad=mean(vals,axis=0)
    print(ad)
    indices_to_keep=[n for n,a in enumerate(ad) if a>1e-6 and a<1.0-1e-6]
    p_new={}
    for k,v in zip(names,vals):
        p_new[k]=array(v)[indices_to_keep]
    return p_new

def add_noise2(p, branch_length):
    sd=sqrt(branch_length)
    p=[a+norm.rvs(scale=sd) for n,a in enumerate(p) ]
    return p

def merge_ps(w, p1,p2):
    #removedprin p1
    #removedprin p2
    res=[w*a+(1-w)*b for a,b in zip(p1,p2)]
    #removedprin res
    return res

def produce_p_matrix(tree, nSNP, clip=True, middle_start=False, allele_dependent=True, fixed_in_outgroup='out'):
    
    ps=uniform.rvs(size=nSNP)
    if middle_start:
        ps=ps*0.2+0.4
    else:
        ps=ps*0.9998+0.0001
    if clip:
        if allele_dependent:
            add_n=add_noise
        else:
            add_n=add_noise3
    else:
        add_n=add_noise2
    ps2=deepcopy(ps)
    p_org=deepcopy(ps)
    res={}
    (child_key1, child_branch1,l1),(child_key2, child_branch2, l2)=find_rooted_nodes(tree)
    
    print(find_rooted_nodes(tree))
    
    if fixed_in_outgroup:
        if child_key1==fixed_in_outgroup:
            tree[child_key2][child_branch2+3]=l1+l2
            tree[child_key1][child_branch1+3]=0
        else:
            tree[child_key2][child_branch2+3]=0
            tree[child_key1][child_branch1+3]=l1+l2
    
    ready_lineages=[(child_key1, child_branch1, ps), (child_key2, child_branch2, ps2)]
    print(ready_lineages)
    started_admixtures={}
    
    while ready_lineages:
        k,b,p = ready_lineages.pop()
        branch_length=get_branch_length(tree,k,b)
        branch=(k,b)
        p=add_n(p,branch_length)
        node=tree[k]
        children=get_children(node)
        if node_is_leaf_node(node):
            res[k]=p
        elif node_is_admixture(node):
            mate=(k,other_branch(b))
            if mate in started_admixtures:
                w=get_admixture_proportion_from_key(tree, k)
                p=merge_ps(w*(1-b)+b*(1-w), p, started_admixtures[mate])
                del started_admixtures[mate]
                k_new=children[0]
                b_new=mother_or_father(tree, k_new, k)
                ready_lineages.append((k_new, b_new, p))
            else:
                started_admixtures[branch]=p
        else:
            for child in children:
                k_new=child
                b_new=mother_or_father(tree, k_new, k)
                ready_lineages.append((k_new,b_new,p))
    return res

def simulate_row_with_binomial(ps, N):
    return binom.rvs(N, ps)/float(N)

def simulate_with_binomial(ps, Ns, p_clip_value=0.01):
    if isinstance(Ns, int):
        Ns={k:Ns for k in ps}
    sims={}
    for k in ps:
        sims[k]=simulate_row_with_binomial(clip(ps[k],p_clip_value, 1.0-p_clip_value), Ns[k])
    return sims

def calculate_covariance_matrix_from_p(ps, nodes=None):
    p_mat=[]
    if nodes is None:
        nodes=list(ps.keys())
    for node in nodes:
        p_mat.append(ps[node])
        
    m=array(p_mat)-mean(p_mat, axis=0)
    return cov(m)
        
    
    