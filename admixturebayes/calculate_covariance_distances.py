from construct_filter_choices import make_filter
from construct_estimator_choices import make_estimator
import numpy as np
from scipy.stats import wishart
from tree_statistics import identifier_file_to_tree_clean, identifier_to_tree_clean
from Rtree_to_covariance_matrix import make_covariance
from reduce_covariance import reduce_covariance, Areduce
from construct_estimator_choices import make_estimator
from copy import deepcopy
import warnings

from tree_to_data import get_xs_and_ns_from_treemix_file, reorder_covariance, order_covariance
import subprocess
import pandas as pd
from Rtree_operations import add_outgroup

def get_average_deviation(scale, df, reps=1000):
    sims=[wishart.rvs(scale=scale/df, df=df) for _ in range(reps)]
    wishs=[wishart.logpdf(cov, df=df, scale=scale/df) for cov in sims]
    dists=[dist(cov, scale) for cov in sims] #ANDREW DEBUG
    return (min(wishs), max(wishs)),(min(dists),max(dists))

def frobenius(a,b):
    return np.linalg.norm(a-b)

def make_wishart_distancer(df):
    def wishart_distance(a,b):
        try:
            ans=(wishart.logpdf(a, df=df, scale=b/df)+wishart.logpdf(b, df=df, scale=a/df))/2
        except np.linalg.LinAlgError as e:
            ans=np.NaN
        return ans
    return wishart_distance

def open_df_file(df_file):
    with open(df_file, 'r') as f:
       df=float(f.readline().rstrip())
    return df

def open_cov_file_admb(cov_mult_file, var_correction, nodes=None):
    res=[]
    with open(cov_mult_file, 'r') as f:
        nodes_from_file=f.readline().split()
        n=len(nodes_from_file)
        for _ in range(n):
            res.append(list(map(float,f.readline().split()[1:])))
        a=f.readline()
        if len(a)>1:
            #removedprin a
            multiplier=float(a.split("=")[1].rstrip())
        else:
            multiplier=1
    res=np.array(res)
    if nodes is None and 's1' in nodes_from_file:
        nodes=['s'+str(i) for i in range(1,len(nodes_from_file)+1)]
    res=reorder_covariance(res, nodes_from_file, nodes)
    if var_correction is not None:
        vc=np.loadtxt(var_correction)
        vc=reorder_covariance(vc, nodes_from_file, nodes)
    else:
        vc=np.zeros(res.shape)
    return res/multiplier, vc, multiplier

def get_true_mat(true_tree_file, nodes):
    scaled_true_tree=identifier_file_to_tree_clean(true_tree_file)
    C=make_covariance(scaled_true_tree, node_keys=nodes)
    return Areduce(C)

def make_A_estimate(xnn, nodes):
    est=make_estimator(reduce_method='average', 
                   variance_correction='unbiased', 
                   indirect_correction=False,
                   nodes=nodes, 
                   arcsin_transform=False, 
                   method_of_weighing_alleles='average_sum', 
                   reducer='',
                   jade_cutoff=1e-5,
                   reduce_also=False,
                   bias_c_weight='default',
                   EM_maxits=0,
                   Indirect_its=0,
                   EM_alpha=0,
                   Indirect_multiplier_s=0,
                   Simulator_fixed_seed=False,
                   initial_Sigma_generator={'default':''},
                   no_repeats_of_cov_est=1,
                   ns=xnn[1], #only necessary if indirect estimation is used.
                   Simulator_from_file='',
                   locus_filter_on_simulated=None,
                   add_variance_correction_to_graph=True,
                   prefix='ccd_',
                   save_variance_correction=True)
    extra_info={}
    cov=est(xnn[0], xnn[1], extra_info=extra_info)
    cov=reorder_covariance(cov, xnn[2], nodes)
    vc=np.loadtxt('ccd_variance_correction.txt')
    vc=reorder_covariance(vc, xnn[2], nodes)
    np.savetxt('ccd_variance_correction.txt', vc)
    var_correction2=np.loadtxt('ccd_variance_correction.txt')
    return cov, var_correction2, extra_info['m_scale']

def load_treemix_matrix(filename, nodes):
    if filename.endswith('.gz'):
        subprocess.call(['gunzip', '-kf', filename])
        filename2='.'.join(filename.split('.')[:-1])
    else:
        filename2=filename
    name_of_file=filename2.split('.')[-1]
    if 'cov' in name_of_file:
        with open(filename2, 'r') as f:
            header=f.readline().split()
            print(header)
            ans=np.zeros((len(header), len(header)))
            #removedprin ans
            for n,lin in enumerate(f.readlines()):
                if len(lin)>1:
                    a=list(map(float, lin.rstrip().split()[1:]))
                    ans[n,:]=a
    a,h =ans, header
    reverse_dic={hi:i for i,hi in enumerate(h)}
    anew=deepcopy(a)
    for n1,k1 in enumerate(nodes):
        for n2,k2 in enumerate(nodes):
            anew[n1,n2]=a[reverse_dic[k1], reverse_dic[k2]]
    return anew

def get_posterior_A_matrices(outfile, add_multiplier=1, nodes=None, outgroup='out', thinning=100):
    a=pd.read_csv(outfile, usecols=['tree','add','layer'])
    b=a.loc[a.layer == 0, :]
    b=b[int(b.shape[0])/2::thinning]
    AmatricesA=[]
    for stree, add in zip(b['tree'], b['add']):
        #removedprin stree
        tree=identifier_to_tree_clean(stree)
        #removedprin pretty_string(tree)
        tree= add_outgroup(tree,  inner_node_name='new_node', to_new_root_length=float(add)*add_multiplier, to_outgroup_length=0, outgroup_name=outgroup)
        cov=make_covariance(tree, node_keys=nodes)
        #removedprin cov
        AmatricesA.append(Areduce(cov))
    return AmatricesA
    

def get_all_dists(nodes_with_outgroup,
                  nodes_without_outgroup,
                  true_tree_file=None, 
                  admB_output_file=None, 
                  admB_covariance=None, 
                  admb_vc=None,
                  snp_data_file=None,
                  treemix_covariance=None,
                  treemix_covariance_without_correction=None,
                  treemix_model_covariances=[],
                  treemix_likelihoods=[],
                  treemix_covse=None,
                  dfs=[],
                  print_dists=True):
    outgroup=(set(nodes_with_outgroup)-set(nodes_without_outgroup)).pop()
    matrices_to_compare={}#contains name:(A/R, print_all/print_min_max_mean, VC True/False, matrix)
    vcs={}
    #multiplier=0.1037936
        
    if admB_covariance is not None:
        RadmbR, RvcR, multiplier = open_cov_file_admb(admB_covariance, admb_vc, nodes_without_outgroup)
        matrices_to_compare['AdmB as input to MCMC']=('R', 'print', False, [RadmbR])
        vcs['R']=RvcR
        if admB_output_file:
            Apost_matsA = get_posterior_A_matrices(admB_output_file, add_multiplier=1.0/multiplier, nodes=nodes_with_outgroup, outgroup=outgroup)
            matrices_to_compare['Posterior samples']=('A','print_min_max_mean', True, Apost_matsA)
    if snp_data_file is not None:
        locus_filter=make_filter('none')
        xnn_tuple=get_xs_and_ns_from_treemix_file(snp_data_file, locus_filter)
        xnn_tuple=order_covariance(xnn_tuple, outgroup=outgroup)
        AadmbA, AvcA, m_scale= make_A_estimate(xnn_tuple, nodes=nodes_with_outgroup)
        matrices_to_compare['AdmB re-estimated in A-space']=('A', 'print', False, [AadmbA])
        vcs['A']=AvcA
    else:
        m_scale=None
    if true_tree_file is not None:
        Atrue_matA = get_true_mat(true_tree_file, nodes_with_outgroup)
        matrices_to_compare['True Matrix']=('A', 'print', True, [Atrue_matA])
        
        
    
    if treemix_covariance is not None:
        AtreemixA=load_treemix_matrix(treemix_covariance, nodes=nodes_with_outgroup)/m_scale
        matrices_to_compare['Treemix input']=('A','print', True, [AtreemixA])
        if snp_data_file is None:
            warnings.warn('because no SNP file is supplied, the treemix matrices cannot be compared to the other', UserWarning)
    if treemix_covariance_without_correction is not None:
        AtreemixVCA= load_treemix_matrix(treemix_covariance_without_correction, nodes=nodes_with_outgroup)/m_scale
        matrices_to_compare['Treemix input without VC']=('A','print', False, [AtreemixVCA])
        if snp_data_file is None:
            warnings.warn('because no SNP file is supplied, the treemix matrices cannot be compared to the other', UserWarning)
    AmodelTreemixAs=[]
    count=0
    for treemix_model_covariance in treemix_model_covariances:
        AmodelTreemixAs.append(load_treemix_matrix(treemix_model_covariance, nodes_with_outgroup)/m_scale)
        if count==0 and snp_data_file is None:
            warnings.warn('because no SNP file is supplied, the treemix matrices cannot be compared to the other', UserWarning)
        count+=1
    if AmodelTreemixAs:
        matrices_to_compare['Treemix fitted matrices']=('A','print', True, AmodelTreemixAs)
    
    if admB_output_file:
        Apost_matsA = get_posterior_A_matrices(admB_output_file, add_multiplier=1.0/multiplier, nodes=nodes_with_outgroup, outgroup=outgroup)
        matrices_to_compare['Posterior samples']=('A','print_min_max_mean', True, Apost_matsA)    
    dist_measures={'Frobenius':frobenius}
    for df in dfs:
        dist_measures['W_'+str(int(df))]=make_wishart_distancer(df)
    
    n_outgroup=next((n for n, e in enumerate(nodes_with_outgroup) if e==outgroup))
    
    if print_dists:
        for space in ['A','R']:
            for dist_name, dist_measure in list(dist_measures.items()):
                for vc in [True,False]:
                    sets_alredy_taken=[]
                    if not (space=='A' and dist_name.startswith('W')) and (dist_name.startswith('W') or vc):
                        print('------ Distance ',dist_name, ' in the distance space ', space + ', vc='+str(vc)+ ' ----------')
                        for matrix_name, (mat_space, summarize_method, vc_corrected, mats) in list(matrices_to_compare.items()):
                            for matrix_name2, (mat_space2, summarize_method2, vc_corrected2, mats2) in list(matrices_to_compare.items()):
                                if (matrix_name,matrix_name2) not in sets_alredy_taken and mat_space<=space and mat_space2<=space and mats and mats2:
                                    mats1,mats2=prepare_mats(mats, mats2, mat_space, mat_space2, vc_corrected, vc_corrected2, space, vc, vcs, n_outgroup)
                                    str_to_print='{:60}'.format(matrix_name+ '><'+matrix_name2)
                                    if summarize_method=='print':
                                        if summarize_method2=='print':
                                            str_to_print+=print_all_print_all(mats1, mats2, dist_measure)
                                        else:
                                            str_to_print+= print_all_print_min_max_mean(mats1, mats2, dist_measure)
                                    elif summarize_method2=='print':
                                        str_to_print+= print_all_print_min_max_mean(mats2, mats1, dist_measure)
                                    print(str_to_print)
                                    sets_alredy_taken.extend([(matrix_name, matrix_name2),(matrix_name2,matrix_name)])
    return matrices_to_compare, vcs


def prepare_mats(mats1, mats2, mat_space1, mat_space2, vc1,vc2, target_space, target_vc, vcs, n_outgroup):
    mats1=[mat+(int(vc1)-int(target_vc))*vcs[mat_space1] for mat in mats1]
    mats2=[mat+(int(vc2)-int(target_vc))*vcs[mat_space2] for mat in mats2]
    if mat_space1=='A' and target_space=='R':
        mats1=[reduce_covariance(mat,n_outgroup) for mat in mats1]
    if mat_space2=='A' and target_space=='R':
        mats2=[reduce_covariance(mat, n_outgroup) for mat in mats2]             
    return mats1, mats2
        
def print_print_all(mat_A, mats_B, dist_m):
    str_to_print=''
    if len(mats_B)>1:
        str_to_print+='\n'
    for mat_B in mats_B:
        str_to_print+='{:6f}'.format(dist_m(mat_A, mat_B))+'\n'
    str_to_print=str_to_print[:-1]
    return str_to_print

def print_all_print_all(mats_A, mats_B, dist_m):
    str_to_print=''
    for mat_A in mats_A:
        str_to_print+=print_print_all(mat_A, mats_B, dist_m)+'\n'
    str_to_print=str_to_print[:-1]
    return str_to_print

def print_print_min_max_mean(mat_A, mats_B, dist_m):
    dists=np.array([dist_m(mat_A, mat_B) for mat_B in mats_B])
    return '{:6f}({:5f},{:5f})'.format(np.mean(dists),np.min(dists),np.max(dists))

def print_all_print_min_max_mean(mats_A, mats_B, dist_m): 
    '''
    A is the print_all and mats_B are the print some
    '''
    str_to_print=''
    if len(mats_A)>1:
        str_to_print+='\n'
    for mat_A in mats_A:
        str_to_print+=print_print_min_max_mean(mat_A, mats_B, dist_m)+'\n'
    str_to_print=str_to_print[:-1]
    return str_to_print

if __name__=='__main__':
    nodes_without_outgroup=['s'+str(i) for i in range(1,10)]
    prefix='../../../../Dropbox/Bioinformatik/AdmixtureBayes/treemix_example3/'
    get_all_dists(nodes_with_outgroup=['out']+nodes_without_outgroup,
                  nodes_without_outgroup=nodes_without_outgroup,
                  true_tree_file=prefix+'_scaled_true_tree.txt', 
                  admB_output_file=prefix+'1_outnew.txt', 
                  admB_covariance=prefix+'_covariance_and_multiplier.txt', 
                  admb_vc=prefix+'_variance_correction.txt',
                  snp_data_file=prefix+'_treemix_in.txt.gz',
                  treemix_covariance=prefix+'new_one2.cov.gz',
                  treemix_covariance_without_correction=None,
                  treemix_model_covariances=[prefix+'new_one.modelcov', prefix+'new_one2.modelcov'],
                  treemix_likelihoods=[],
                  treemix_covse=None,
                  dfs=[1000,10000])
    
        
    