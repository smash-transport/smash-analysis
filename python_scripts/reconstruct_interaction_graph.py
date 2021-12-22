#!/usr/bin/python3

"""
reconstruct_interaction_graph.py

This script contains all the functions necessary to reconstruct the
interaction graph of the particles, the maximum value of the number
of mother particles for each particle and its offspring, and the final
states of the resonance.
"""

from collections import defaultdict, Counter


class Particle:
    """A particle node of the interaction graph."""
    def __init__(self):
        self.mothers = []
        self.daughters = []
        self.pdgid = None
        self.pos = []
        self.mom = []

    def __repr__(self):
        return 'Particle(pdgid={}, mothers={}, daughters={})'.format(self.pdgid,
                 self.mothers, self.daughters)


def indegrees(nodecay_list, graph):
    """Calculate the number of mothers for each particle."""
    c = Counter()
    for k in list(graph.keys()):
        if set(graph[k].daughters) - set(graph[k].mothers) == set([]): # list the particles that do not decay
           nodecay_list.append(k)
        for d in graph[k].daughters:
            c[d] += 1
    return c

def max_indegrees_bottom_up(final_list, graph):
    """Calculate the maximum number of mothers going upwards from the given
    particles."""
    nodecay_list = []
    indeg = indegrees(nodecay_list, graph)
    starts = set(nodecay_list).union(set(final_list)) # starts from the particles either in the final state or removed in the thermal bubble
    removed = set(nodecay_list) - set(final_list) # list the particles removed in the thermal bubble
    stack = [(s, indeg[s]) for s in starts]
    maxindeg_bookkeeper = {}
    while stack:
        node, max_indeg = stack.pop()
        if node in removed:
           max_indeg = float('inf') # If a particle is removed in the thermal bubble, its mother, grand mother, grand grand mother, etc. shouldn't be observed
        old_maxindeg = maxindeg_bookkeeper.get(node, 0)
        max_indeg = max(max_indeg, old_maxindeg)
        maxindeg_bookkeeper[node] = max_indeg
        mothers = set(graph[node].mothers)
        daughters = set(graph[node].daughters)
        # An elastic interaction between A and B adds three cycles to the graph:
        # A -> A, B -> B and A -> B -> A.
        # To break those cycles, we ignore any mothers that are also our daughters.
        # This avoids an infinite loop while still visiting all particles.
        for mother in mothers - daughters:
            max_indeg = max(max_indeg, indeg[mother])
            stack.append((mother, max_indeg))
    return maxindeg_bookkeeper

def count_interaction(interaction_graph, datablock):
    """Count the processes in all interactions.

    Fills the `interaction_graph`."""
    assert datablock['type'] == 'i'
    final_ids = datablock['outgoing']['id']
    for p in datablock['incoming']:
        mother = interaction_graph[p['id']]
        #assert mother.daughters == [], "{}".format(mother.daughters)
        mother.daughters.extend(final_ids)
        mother.pdgid = p['pdgid']
        mother.pos = p['r']
        mother.mom = p['p']
        for fp in datablock['outgoing']:
            daughter = interaction_graph[fp['id']]
            daughter.mothers.append(p['id'])
            daughter.pdgid = fp['pdgid']
            daughter.pos = fp['r']
            daughter.mom = fp['p']

def get_daughters(interaction_graph, pid):
    """Take the given particle ID and determine the corresponding final state."""
    stack = interaction_graph[pid].daughters[:]
    daughters = []
    while stack:
        current = stack.pop()
        current_daughters = set(interaction_graph[current].daughters)
        # Ignore elastic interactions, they don't change anything.
        # (This also breaks all cycles in the interaction graph.)
        current_daughters -= set(interaction_graph[current].mothers)
        if current_daughters:
            stack.extend(current_daughters)
        else:
            # The particle is stable and will be there
            # after the full decay of the mother.
            daughters.append(current)
    return daughters
