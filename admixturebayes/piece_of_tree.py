class piece(object):
    
    def __init__(self, start_lattitude, end_lattitude, start_distance, end_distance, child_key, child_branch, parent_key):
        self.start_lattitude=start_lattitude
        self.end_lattitude=end_lattitude
        self.start_distance=start_distance
        self.end_distance=end_distance
        self.child_key=child_key
        self.child_branch=child_branch
        self.parent_key=parent_key
        if end_lattitude is not None and start_lattitude>end_lattitude:
            self.direction='to_leaves'
        else:
            self.direction='to_root'
        
    def __str__(self):
        return ', '.join(map(str,[(self.child_key,self.child_branch, self.parent_key), self.start_lattitude, self.end_lattitude, self.start_distance, self.end_distance, self.direction]))
    
    def get_branch_key(self):
        return (self.child_key, self.child_branch)
    
    def get_coverage(self):
        return self.start_distance, self.end_distance
    
    def contains_distance(self, distance):
        if self.end_distance is None:
            return distance>=self.start_distance
        return distance<=self.end_distance and distance>=self.start_distance
    
    def within_distance(self, distance):
        return self.start_distance<=distance
    
    def get_start_distance(self):
        return self.start_distance
    
    def get_leaf_and_root_sided_length(self, distance):
        if self.end_distance is None:
            return distance-self.start_distance, None
        if self.direction=='to_leaves':
            return self.end_distance-distance, distance-self.start_distance
        else:
            return distance-self.start_distance, self.end_distance-distance
        
    def get_lattitude(self, distance):
        if self.end_distance is None:
            return self.start_lattitude+distance
        elif self.direction == 'to_leaves':
            return self.start_lattitude-distance
        else:
            return self.start_lattitude+distance