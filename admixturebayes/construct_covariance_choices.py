from tree_to_data import (emp_cov_to_file,
                          get_xs_and_ns_from_treemix_file, order_covariance, reorder_reduced_covariance)
from copy import deepcopy
from numpy import loadtxt, savetxt
from construct_estimator_choices import make_estimator
def rescale_empirical_covariance(m):    
    n=m.shape[0]
    actual_trace=m.trace()
    max_expected_trace=n*(n+1)/2-1
    multiplier= max_expected_trace/actual_trace
    
    return m*multiplier, multiplier

def empirical_covariance_wrapper_directly(snp_data_file, **kwargs):
    xnn_tuple=get_xs_and_ns_from_treemix_file(snp_data_file)
    return xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs)

def xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs):
    est_args=kwargs['est']
    xnn_tuple=order_covariance(xnn_tuple, outgroup=est_args['reducer'])
    xs,ns,names=xnn_tuple

    est= make_estimator(reduce_method='outgroup', **est_args)
    extra_info_dic={}
    cov=est(xs,ns, extra_info_dic)
    cov=reorder_reduced_covariance(cov, names, est_args['nodes'], outgroup=est_args['reducer'])
    if ('add_variance_correction_to_graph' in est_args and
        est_args['add_variance_correction_to_graph'] and
        'save_variance_correction' in est_args and
        est_args['save_variance_correction']):
        filename=est_args['prefix']+'variance_correction.txt'
        vc=loadtxt(filename)
        vc=reorder_reduced_covariance(vc, names, est_args['nodes'], outgroup=est_args['reducer'])
        savetxt(filename, vc)
    if 'm_scale' in extra_info_dic:
        if 'import' in kwargs:
            with open(kwargs['import'], 'w') as f:
                txt=str(extra_info_dic['m_scale'])
                assert '\0' not in txt, 'binary content in m_scale file'
                f.write(txt)
        if 'return_also_mscale' in kwargs and kwargs['return_also_mscale']:
            return cov, extra_info_dic['m_scale']
    return cov

def normaliser_wrapper(covariance, **kwargs):
    return rescale_empirical_covariance(covariance)

dictionary_of_transformations={
    (6,8):empirical_covariance_wrapper_directly,
    (8,9):normaliser_wrapper
    }

dictionary_of_reasonable_names={
    1:'number_of_leaves',
    2:'leaves_admixtures',
    3:'true_tree',
    4:'true_tree_with_outgroup',
    5:'scaled_true_tree',
    6:'SNP_data',
    7:'covariance',
    8:'covariance_without_reduce_name',
    9:'covariance_and_multiplier'}

def write_one_line_to_file(filename, value):
    with open(filename,'w') as f:
        f.write(value)

def write_two_lines_to_file(filename, value1, value2):
    with open(filename, 'w') as f:
        f.write(value1+'\n'+value2)

def save_stage(value, stage_number, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes, filename=None):
    if filename is None:
        save_word=dictionary_of_reasonable_names[stage_number]
        filename=prefix+save_word+'.txt'
    if stage_number==8:
        emp_cov_to_file(value, filename, after_reduce_nodes)
    else:
        emp_cov_to_file(value[0], filename, after_reduce_nodes)
        with open(filename, 'a') as f:
            f.write('multiplier='+str(value[1]))


def get_covariance(input, full_nodes=None,
                   p=0.5,
                   reduce_covariance_node=None,
                   blocksize_empirical_covariance=100,
                   save_stages=list(range(1,6))+list(range(7,10)),
                   prefix='tmp',
                   estimator_arguments={}):

    kwargs={}
    kwargs['p']=p
    kwargs['full_nodes']=full_nodes
    before_added_outgroup_nodes=deepcopy(full_nodes)
    after_reduce_nodes=deepcopy(full_nodes)
    after_reduce_nodes.remove(reduce_covariance_node)
    print(before_added_outgroup_nodes)
    print(after_reduce_nodes)

    print(reduce_covariance_node)

    kwargs['reduce_covariance_node']=reduce_covariance_node
    kwargs['after_reduce_nodes']=after_reduce_nodes
    kwargs['before_added_outgroup_nodes']=before_added_outgroup_nodes
    kwargs['blocksize_empirical_covariance']=blocksize_empirical_covariance
    kwargs['pks']={}
    kwargs['add_file']=prefix+'true_add.txt'
    kwargs['import']=prefix+'m_scale.txt'
    kwargs['est']=estimator_arguments

    stages_to_go_through = [6,8,9]
    #makes a necessary transformation of the input(if the input is a filename or something).
    statistic = input
    for stage_from, stage_to in zip(stages_to_go_through[:-1], stages_to_go_through[1:]):
        transformer_function=dictionary_of_transformations[(stage_from, stage_to)]
        statistic=transformer_function(statistic, **kwargs)
        if stage_to in save_stages:
            save_stage(statistic, stage_to, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)

    return statistic

def read_one_line(filename):
    with open(filename, 'r') as f:
        res=f.readline().rstrip()
    return res
