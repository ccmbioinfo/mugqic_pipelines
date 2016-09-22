#!/usr/bin/env python

"""
Visualize job dependencies

See dependency_graph below
"""

import networkx


###################################################################################################
# CREATE GRAPH
###################################################################################################

def lookup_id(job_logs, id):
    """
    Find the JobLog with this id
    """
    return next(log for log in job_logs if log.job_id == id)


def create_dependency_graph(job_logs, dependency_first=True):
    """
    Read the job_logs and return a dependency graph, where each node is a JobLog

    :param job_logs: list of JobLogs
    :param dependency_first: whether edges go from dependency to depender
    :return: networkx.DiGraph
    """
    graph = networkx.DiGraph()

    for log in job_logs:
        for id in log.job_dependencies:
            # Look the log associated with this id
            dependency = lookup_id(job_logs, id)

            if dependency_first:
                graph.add_edge(dependency, log)
            else:
                graph.add_edge(log, dependency)

    return graph


###################################################################################################
# RENDER GRAPH
###################################################################################################

def recursive_subgraph(nodes, graph):
    return set().union(*({node} | networkx.descendants(graph, node) for node in nodes))


def sort_by_height(nodes, graph):
    """
    Sort the nodes in graph by the longest distance to a leaf

    Take the nodes into a subgraph, and call 'topological_sort' on that new graph
    :param nodes: list of JobLogs
    :param graph: networkx.Digraph
    :return: sorted list of JobLogs
    """
    subgraph = graph.subgraph(recursive_subgraph(nodes, graph))
    return [node for node in reversed(networkx.topological_sort(subgraph)) if node in nodes]


def render_node(log, graph, text, covered_ids, num_indents=0):
    """
    :param log: JobLog
    :param graph: networkx.Digraph
    :param text: text to append to
    :param covered_ids: set of job ids
    :param num_indents: number of tabs to indent this node
    :return: text, set of JobLogs that were covered in this render
    """
    dependency_string = '(' + str(log.job_dependencies) + ')' if len(log.job_dependencies) > 0 else ''
    #TODO: status string
    # status_string = ' <- Running' if hasattr(log, 'status') and log.status == 'INACTIVE' else ''

    text += log.job_id + '|  ' * (num_indents + 1) + log.job_name + ' ' + dependency_string + '\n'
    # text += log.job_id + '|\t' * (num_indents + 1) + log.job_name + ' ' + dependency_string + status_string + '\n'

    covered_ids.add(log.job_id)

    # for successor in graph.successors(log):
    for successor in sort_by_height(graph.successors(log), graph):
        if set(successor.job_dependencies) <= covered_ids:
            # Render this node too
            text, covered_ids = render_node(successor, graph, text, covered_ids, num_indents + 1)

    return text, covered_ids


def dep_graph_to_text(graph):
    """
    :param graph: networkx.Digraph
    :return: str
    """
    text = ''
    roots = [node for node in graph.nodes() if len(graph.predecessors(node)) == 0]
    covered_ids = set()

    for index, root in enumerate(sort_by_height(roots, graph)):
    # for index, root in enumerate(roots):
        text, covered_ids = render_node(root, graph, text, covered_ids=covered_ids, num_indents=0)

    return text


###################################################################################################
# MAIN
###################################################################################################

def dependency_graph(job_logs, dependency_first=True):
    """
    :param dependency_first: whether the leftmost jobs should be the dependencies of other jobs
    :return: text representing dependency graph
    """
    graph = create_dependency_graph(job_logs, dependency_first)
    return dep_graph_to_text(graph)
