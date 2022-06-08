from copy import deepcopy
from Rtree_operations import get_trivial_nodes
from tree_to_data import unzip

def get_nodes(arguments, input_file, outgroup_name, reduce_node, backup_number=8):
    ''' The outgroup_name is only used for simulation purposes and reduce_node is the important one
    that is used when analysing the admixture graphs.  '''
    if not arguments[0]:#this means that we should use the input file for nodes
        if ';' in input_file:
            nodes=get_trivial_nodes(len(input_file.split('-')[0].split('.')))
        elif '.' in input_file:
            nodes=read_one_line(input_file)
        elif ',' in input_file:
            nodes=get_trivial_nodes(int(input_file[1:].split(',')[0]))
        else:
            nodes=get_trivial_nodes(int(input_file))
    else:
        nodes=arguments
    before_added_outgroup=deepcopy(nodes)
    reduced_nodes=deepcopy(nodes)
    if outgroup_name in nodes:
        before_added_outgroup.remove(outgroup_name)
    elif outgroup_name:
        nodes.append(outgroup_name)
    if reduce_node in reduced_nodes:      
        reduced_nodes.remove(reduce_node)
    if reduce_node and reduce_node not in nodes:
        nodes.append(reduce_node)
    return before_added_outgroup, nodes, reduced_nodes
        
def read_one_line(filename):
    if filename.endswith('.gz'):
        filename=unzip(filename)
    with open(filename, 'r') as f:
        return f.readline().rstrip().split()
        