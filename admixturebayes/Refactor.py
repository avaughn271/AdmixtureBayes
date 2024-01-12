import random

class Node(object):
    def __init__(self):
        self.left_parent = #left parent, this is a node. If this is none, then this nodeis the root.
        self.right_parent = #right parent, this is a node or None if it is a divergence node/root
        self.left_child = #left_child, this is none if this is a leaf node.
        self.right_child = #right_child, this is none if this is a leaf node/admixture node.
        self.left_branch = #number, branch length of branch to left parent. or none if this node is the root.
        self.right_branch = #number, branch length of branch to right parent or None if this node is divergence node/root.
        self.left_admix_prop = #proportion that came from the left parent_branch or None if not admixture node.

class AdmixtureGraph(object):
    def __init__(self):
        self.nodes = {} # list of all nodes in the ARG is an index:node object pair.
        self.numberofadmixes = -1
        self.outgroupdistance = -1

#seven proposals

def rescale_alladmixtures(graph): # have not checked.
	for i in graph.nodes:
		nodei = graph.nodes[i]
		if isadmixnode(nodei):
			graph.nodes[nodei].left_admix_prop = random.uniform(0, 1)
	return graph

def rescale_outgroup(graph, stdev): # have not checked.
	graph.outgroupdistance += random.gauss(0, stdev)

def branch_random_walk(graph, stdev): # have not checked.
	for i in graph.nodes:
		nodei = graph.nodes[i]
		if nodei.left_branch != None:
			graph.nodes[nodei].left_branch += random.uniform(0, stdev)
		if nodei.right_branch != None:
			graph.nodes[nodei].right_branch += random.uniform(0, stdev)
	return graph

#1[A]-add_ad
#2[A]-remove_ad
#3_hard_part-node-sliding
#4 done-random walk on branch lengths (maybe change to just one branch)
#5[A]done-resample admixture proportions
#6 done-randomwalk on branch to outgroup
#7[don't do]-randomwalk inside the nullspace

def OneStep(graph):
    if ind == 0:
		print("PROBLME")
	elif ind == 1:
		if len(self.arg.rec) != 0:
			self.remove_recombination()
	elif ind == 2:
		self.add_recombination()
	elif ind == 3:
		self.adjust_times()
	elif ind == 4:
		self.adjust_breakpoint()
	elif ind == 5:
		self.kuhner()
	elif ind == 6:
		print("PROBLME44")
	self.detail_acceptance[ind][0] += 1