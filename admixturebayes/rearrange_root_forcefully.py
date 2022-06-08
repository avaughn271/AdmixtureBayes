from Rtree_operations import get_all_admixture_meetings, rearrange_root, remove_admix
from tree_statistics import get_timing

def rearrange_root_force(tree, new_outgroup):
    '''
    Like rearrange_root this changes the outgroup by reversing branches. However, before doing so, all obstacles are removed.
    '''
    if tree[new_outgroup][0]=='r' or tree[new_outgroup][1]=='r':
        #The root was already in the requested location so no rearranging performed
        return tree
    relevant_admixture_keys=get_all_admixture_meetings(tree, new_outgroup)
    key_timings=get_timing(tree)
    timed_admixtures=[(key, key_timings[key]) for key in relevant_admixture_keys]
    s_adm=sorted(timed_admixtures, key=lambda  x: x[1], reverse=True)
    for admixture_to_remove,_ in s_adm:
        tree=remove_admix(tree, admixture_to_remove, 1)[0]
    return rearrange_root(tree, new_outgroup)
