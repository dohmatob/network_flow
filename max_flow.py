"""
Home-made algorithms for single-commodity max-flow in directed networks.
"""

#Author: Dohmatob Elvis <gmdopp@gmail.com>

import networkx as nx
import numpy as np
from scipy import optimize


class FlowNetWork(nx.DiGraph):
    def __init__(self, formulation="lp", *args, **kwargs):
        super(FlowNetWork, self).__init__(self, *args, **kwargs)
        self.formulation = formulation
        self.flow = {}

    def add_edge(self, src, dst, capacity=0.):
        edge = src, dst
        redge = dst, src
        super(FlowNetWork, self).add_edge(src, dst, capacity=capacity)
        self.flow[edge] = 0.
        if self.formulation == "ff":
            super(FlowNetWork, self).add_edge(dst, src, capacity=0.)
            self.flow[redge] = 0.

    def _path_contains_edge(self, path, edge):
        """XXX Very very inefficient method!"""
        for e, _ in path:
            if e == edge:
                return True
        return False

    def find_augmenting_path(self, s, t, path):
        """Returns augmenting path, or None if doesn't exist.

        Parameters
        ----------
        s: node
           Source node.

        t: node
            Sink node.
        """
        if s == t:
            return path
        for sink, data in self.edge[s].items():
            edge = s, sink
            capacity = data["capacity"]
            residual = capacity - self.flow[edge]
            if residual > 0. and not self._path_contains_edge(path, edge):
                result = self.find_augmenting_path(sink, t,
                                                   path + [(edge, capacity)])
                if result is not None:
                    return result

    def _reversed_edge(self, edge):
        return edge[::-1]

    def max_flow(self, s, t, max_iter=-1):
        """Finds Max flow between source 's' and sink 't' in a directed network.

        Parameters
        ----------
        s: node
           Source node.

        t: node
            Sink node.

        max_iter: int, optional (default None)
           Maximum number of iterations after which algorithm will be abrupted.
           A negative value means "run the algorithm until convergence".

        Notes
        -----
        Note that thes algorithm is not guaranteed to terminate.
        """
        if self.formulation == "lp":
            # precompute node and edge indices for subsequent referencing
            node_keys = {}
            for n, node in enumerate(self.nodes()):
                node_keys[node] = n
            edge_keys = {}
            for e, edge in enumerate(self.edges()):
                edge_keys[edge] = e

            # Build constraints. The problem to solve is:
            # Maximize sum {x(s, u): (s, u) in E}
            # subject to:
            #            conservation law: in_flow(v) = out_flow(v)
            #                          for all internal nodes v.
            #            capacity constraints: 0 <= x(e) <= c(e)
            #                          for all edges e.
            A = np.zeros((G.number_of_nodes(), G.number_of_edges()))
            c = np.zeros(A.shape[1])
            bounds = np.zeros((len(c), 2))
            for n, node in enumerate(self.node):
                if node == t:
                    continue
                for out_edge in G.out_edges(node):
                    e = edge_keys[out_edge]
                    bounds[e][1] = G.get_edge_data(*out_edge)["capacity"]
                    if node == s:
                        c[e] = 1.
                    else:
                        A[n, e] = -1.
                if node != s:
                    for in_edge in G.in_edges(node):
                        A[n, edge_keys[in_edge]] = 1.
            b = np.zeros_like(A[:, 0])

            # compute flow by solving LCP
            out = optimize.linprog(-c, A_eq=A, b_eq=b, bounds=bounds)
            print out
            print "=" * 80
            flow = -out.fun
            for edge, e in edge_keys.items():
                self.flow[edge] = out.x[e]
        elif self.formulation == "ff":
            # Ford-Fulkerson algorithm adapted from wikipedia page:
            # https://en.wikipedia.org/wiki/Ford%E2%80%93Fulkerson_algorithm
            # Note that this algorithm is not guaranteed to terminate.
            k = max_iter
            path = self.find_augmenting_path(s, t, [])
            while path is not None and k != 0:
                k -= 1
                residuals = [capacity - self.flow[edge]
                             for edge, capacity in path]
                min_flow = min(residuals)
                for edge, _ in path:
                    redge = self._reversed_edge(edge)
                    self.flow[edge] += min_flow
                    self.flow[redge] -= min_flow
                path = self.find_augmenting_path(s, t, [])
            flow = sum(self.flow[(s, sink)]
                       for sink, _ in self.edge[s].items())

        print "Max flow: %g." % flow
        print "Realization:"
        for edge in G.edges_iter():
            src, dst = edge
            amnt = G.flow[edge]
            if amnt != 0.:
                if amnt < 0.:
                    amnt *= -1.
                    src, dst = dst, src
                print "\tSend %g units from %s to %s." % (amnt, src, dst)
        return flow

if __name__ == "__main__":
    # build the network
    G = FlowNetWork(formulation="ff")
    which = "CLRS"
    if which == "CLRS":
        G.add_edge("s", "1", capacity=16)
        G.add_edge("s", "2", capacity=13)
        G.add_edge("1", "2", capacity=10)
        G.add_edge("2", "1", capacity=4)
        G.add_edge("1", "3", capacity=12)
        G.add_edge("3", "2", capacity=9)
        G.add_edge("2", "4", capacity=14)
        G.add_edge("4", "3", capacity=7)
        G.add_edge("3", "t", capacity=20)
        G.add_edge("4", "t", capacity=4)
    elif which == "Orlin":
        G.add_edge("s", "1", capacity=10)
        G.add_edge("1", "t", capacity=8)
        G.add_edge("s", "2", capacity=6)
        G.add_edge("2", "t", capacity=10)
        G.add_edge("1", "2", capacity=1)
    elif which == "wikipedia":
        G.add_edge('s', 'o', capacity=3)
        G.add_edge('s', 'p', capacity=3)
        G.add_edge('o', 'p', capacity=2)
        G.add_edge('o', 'q', capacity=3)
        G.add_edge('p', 'r', capacity=2)
        G.add_edge('r', 't', capacity=3)
        G.add_edge('q', 'r', capacity=4)
        G.add_edge('q', 't', capacity=2)
    else:
        raise ValueError("Unknown problem: %s" % which)

    # plotting
    import matplotlib.pyplot as plt
    edge_labels = nx.get_edge_attributes(G, 'capacity')
    pos = nx.spring_layout(G)
    nx.draw_networkx(G, pos)
    nx.draw_networkx_edge_labels(G, pos, labels=edge_labels)
    plt.show()

    # compute the max-flow
    G.max_flow('s', 't')
