from Rtree_operations import make_consistency_checks

def check(tree, pks={}):
    result, fails = make_consistency_checks(tree)
    if not result:
        stringed_message=''
        print(fails)
        for test_key, (test_result, test_message) in list(fails.items()):
            if test_result:
                stringed_message+=test_key+": TRUE"+'\n'
            else:
                stringed_message+=test_key+": FALSE"+'\n'
                stringed_message+="    "+test_message+'\n'
        print(stringed_message)
        assert result, "The tree was inconsistent :" + '\n'+stringed_message