import nblist
import numpy as np
import time
import freud

nblist.set_num_threads(4)
print(nblist.get_max_threads())

num_particle = 3000
domain_size = [50.0, 50.0, 50.0]
angle = [90.0, 90.0, 90.0]
shape = (num_particle, 3)
inputs = np.random.random(shape) * domain_size
cutoff = 6.0

box = nblist.Box(domain_size, angle)
nb = nblist.NeighborList("Octree-CPU")    # Linked_Cell-GPU  Linked_Cell-CPU  Octree-GPU  Octree-CPU  Hash-GPU  Hash-CPU
nb.build(box, inputs, cutoff)
pairs = nb.get_neighbor_pair()
lists = nb.get_neighbor_list()
print(pairs[0:9])
print(lists[0:3])
