from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from tree_to_data import (file_to_emp_cov,
                         emp_cov_to_file,
                          get_xs_and_ns_from_treemix_file, order_covariance, reorder_reduced_covariance)
from copy import deepcopy
from numpy import loadtxt, savetxt
from construct_estimator_choices import make_estimator
from math import log
def rescale_empirical_covariance(m, normalizer= ['min', 'max']):    
    if not isinstance(normalizer, str):
        normalizer=normalizer[0]
    
    n=m.shape[0]
    actual_trace=m.trace()
    min_expected_trace=log(n)/log(2)*n
    max_expected_trace=n*(n+1)/2-1
    
    if normalizer=='min':
        multiplier= min_expected_trace/actual_trace
    elif normalizer=='max':
        multiplier= max_expected_trace/actual_trace
    else:
        assert False, 'normalizer not set properly. normalizer= '+str(normalizer)
    
    return m*multiplier, multiplier

def empirical_covariance_wrapper_directly(snp_data_file, **kwargs):
    xnn_tuple=get_xs_and_ns_from_treemix_file(snp_data_file, kwargs['locus_filter'])
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
    return rescale_empirical_covariance(covariance, normalizer=kwargs['scale_goal'])

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
                   outgroup_name=None,
                   reduce_covariance_node=None,
                   theta=0.4, sites=500000,
                   treemix_file=None,
                   blocksize_empirical_covariance=100,
                   save_stages=list(range(1,6))+list(range(7,10)),
                   prefix='tmp',
                   final_pop_size=100.0,
                   filter_on_outgroup=False,
                   locus_filter=None,
                   estimator_arguments={},
                   verbose_level='normal'):

    if treemix_file is None:
        treemix_file=prefix+'treemix'

    kwargs={}
    kwargs['p']=p
    kwargs['outgroup_name']=outgroup_name
    kwargs['full_nodes']=full_nodes
    before_added_outgroup_nodes=deepcopy(full_nodes)
    after_reduce_nodes=deepcopy(full_nodes)
    if outgroup_name is not None and outgroup_name in before_added_outgroup_nodes:
        before_added_outgroup_nodes.remove(outgroup_name)
    if reduce_covariance_node is not None and reduce_covariance_node in after_reduce_nodes:
        after_reduce_nodes.remove(reduce_covariance_node)
    kwargs['reduce_covariance_node']=reduce_covariance_node
    kwargs['after_reduce_nodes']=after_reduce_nodes
    kwargs['before_added_outgroup_nodes']=before_added_outgroup_nodes
    kwargs['theta']=theta
    kwargs['sites']=sites
    kwargs['treemix_file']=treemix_file
    kwargs['blocksize_empirical_covariance']=blocksize_empirical_covariance
    kwargs['pks']={}
    kwargs['final_pop_size']=final_pop_size
    kwargs['add_file']=prefix+'true_add.txt'
    kwargs['import']=prefix+'m_scale.txt'
    kwargs['scale_goal']='max'
    kwargs['filter_on_outgroup']=filter_on_outgroup
    kwargs['est']=estimator_arguments
    kwargs['locus_filter']=locus_filter

    stages_to_go_through = [6,8,9]
    #makes a necessary transformation of the input(if the input is a filename or something).
    statistic = input
    for stage_from, stage_to in zip(stages_to_go_through[:-1], stages_to_go_through[1:]):
        transformer_function=dictionary_of_transformations[(stage_from, stage_to)]
        statistic=transformer_function(statistic, **kwargs)
        if stage_to in save_stages:
            save_stage(statistic, stage_to, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)

    return statistic

def read_multiplier(input):
    with open(input, 'r') as f:
        last_line=f.readlines()[-1]
        return float(last_line.split("=")[1])

def read_covariance_matrix(input, nodes):
    if isinstance(input, str):
        return file_to_emp_cov(input, nodes=nodes)
    else:
        return input

def read_tree(input, nodes):
    if isinstance(input, str):
        if not ';' in input:
            input=read_one_line_skip(filename=input)
            return identifier_to_tree_clean(input, leaves=generate_predefined_list_string(deepcopy(nodes)))
        else:
            return identifier_to_tree_clean(input, leaves=generate_predefined_list_string(deepcopy(nodes)))
    else:
        return input

def read_one_line_skip(filename):
    with open(filename, 'r') as f:
        lines=f.readlines()
        if len(lines[-1])>3:
            return lines[-1].rstrip()
        else:
            return lines[-2].rstrip()

def read_one_line(filename):
    with open(filename, 'r') as f:
        res=f.readline().rstrip()
    return res
