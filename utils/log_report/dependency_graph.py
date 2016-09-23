#!/usr/bin/env python

"""
visualize job dependencies

see dependency_graph below
"""

import networkx


###################################################################################################
# create graph
###################################################################################################

def lookup_id(job_logs, id):
    """
    find the joblog with this id
    """
    return next(log for log in job_logs if log.job_id == id)


def create_dependency_graph(job_logs, dependency_first=True):
    """
    read the job_logs and return a dependency graph, where each node is a joblog

    :param job_logs: list of joblogs
    :param dependency_first: whether edges go from dependency to depender
    :return: networkx.digraph
    """
    graph = networkx.DiGraph()

    for log in job_logs:
        graph.add_node(log)

        for id in log.job_dependencies:
            # look the log associated with this id
            dependency = lookup_id(job_logs, id)

            if dependency_first:
                graph.add_edge(dependency, log)
            else:
                graph.add_edge(log, dependency)

    return graph


###################################################################################################
# render graph
###################################################################################################

def height(node, graph):
    try:
        return max([1 + height(succ, graph) for succ in graph.successors(node)])
    except ValueError:
        return 0


def sort_by_height(nodes, graph):
    return sorted(nodes, key=lambda n: height(n, graph))


def render_dependencies(node):
    return ' (' + str(node.job_dependencies) + ')' if len(node.job_dependencies) > 0 else ''


def render_status(node):
    status_to_display = dict({
        'SUCCESS' : ' + COMPLETED',
        'FAILED' : ' - FAILED',
        'ACTIVE' : ' <- Active',
        'INACTIVE' : ' -> Pending'

    })
    return status_to_display[node.status]


def render_node(node, num_indents):
    """
    Render only this node with the given indentation

    :param node: JobLog
    :param num_indents: between job id and job name
    :return: string representing this node
    """
    return node.job_id +\
           num_indents * '|  ' +\
           node.job_name +\
           render_dependencies(node) +\
           render_status(node) + '\n'


def render_nodes(nodes, graph, rendered_nodes=set(), num_indents=1):
    """
    Render then nodes in order of smallest associated text to largest

    :param nodes:
    :param graph:
    :param rendered_nodes: set of nodes that have been displayed already
    :param num_indents: number of indentations between <job_id> and <job_name>
    :return: str representing nodes
    """
    node_texts = []

    for node in sort_by_height(nodes, graph):
    # for node in nodes:
        # We will only render the node if all of its predecessors have been rendered above it
        if set(graph.predecessors(node)) <= rendered_nodes:
            rendered_nodes.add(node)

            node_texts.append(render_node(node, num_indents) +
                              render_nodes(graph.successors(node), graph, rendered_nodes, num_indents + 1))

    return ''.join(sorted(node_texts, key=lambda t: len(t.split())))


def dep_graph_to_text(graph):
    """
    :param graph: networkx.digraph
    :return: str
    """
    text = ''
    roots = [node for node in graph.nodes() if len(graph.predecessors(node)) == 0]
    covered_logs = set()

    return render_nodes(roots, graph)


###################################################################################################
# main
###################################################################################################

def dependency_graph(job_logs, dependency_first=True):
    """
    :param dependency_first: whether the leftmost jobs should be the dependencies of other jobs
    :return: text representing dependency graph
    """
    graph = create_dependency_graph(job_logs, dependency_first)
    return dep_graph_to_text(graph)
