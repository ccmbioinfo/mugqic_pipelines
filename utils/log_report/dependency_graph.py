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


def create_dependency_graph(job_logs, dependencies_first=True):
    """
    read the job_logs and return a dependency graph, where each node is a joblog

    :param job_logs: list of joblogs
    :param dependencies_first: whether edges go from dependency to depender
    :return: networkx.digraph
    """
    graph = networkx.DiGraph()

    for log in job_logs:
        graph.add_node(log)

        for id in log.job_dependencies:
            # look the log associated with this id
            dependency = lookup_id(job_logs, id)

            if dependencies_first:
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
        'SUCCESS' : '',
        'FAILED' : ' - FAILED',
        'ACTIVE' : ' <- ACTIVE',
        'INACTIVE' : ''

    })
    return status_to_display[node.status] if hasattr(node, 'status') else ''


def render_node(node, num_indents):
    """
    Render only this node with the given indentation

    Resulting string is of the following form:
    <job_id>    <job_name> [dependencies] [status]

    :param node: JobLog
    :param num_indents: between job id and job name
    :return: string representing this node
    """
    return node.job_id +\
           num_indents * '|  ' +\
           node.job_name +\
           render_dependencies(node) +\
           render_status(node) + '\n'


def render_nodes(nodes, graph, rendered_nodes=None, num_indents=1):
    """
    Render nodes as an indented list

    First, order the nodes by height and then render each one.
    When assembling the final string, display the shortest texts first

    :param nodes: list of JobLogs
    :param graph: networkx.DiGraph
    :param rendered_nodes: set of nodes that have been displayed already
    :param num_indents: number of indentations between <job_id> and <job_name>
    :return: str representing nodes
    """
    # If no rendered_nodes are given, then assume the empty set
    rendered_nodes = set() if rendered_nodes == None else rendered_nodes

    # List of strings associated with each node
    node_strings = []

    for node in sort_by_height(nodes, graph):
        # We will only render the node if all of its predecessors have been rendered above it
        if set(graph.predecessors(node)) <= rendered_nodes:
            rendered_nodes.add(node)

            # Render this node and as many of its successors as possible
            node_strings.append(render_node(node, num_indents) +
                              render_nodes(graph.successors(node), graph, rendered_nodes, num_indents + 1))

    # Make sure that the shortest strings are listed first
    return ''.join(sorted(node_strings, key=lambda s: len(s.split())))


def dep_graph_to_text(graph):
    """
    :param graph: networkx.DiGraph
    :return: str
    """
    text = ''
    roots = [node for node in graph.nodes() if len(graph.predecessors(node)) == 0]
    covered_logs = set()

    return render_nodes(roots, graph)


###################################################################################################
# main
###################################################################################################

def dependency_graph(job_logs, dependencies_first=True):
    """
    :param dependencies_first: whether the leftmost jobs should be the dependencies of other jobs
    :return: text representing dependency graph
    """
    graph = create_dependency_graph(job_logs, dependencies_first)

    if dependencies_first:
        header = '# Dependency jobs listed first:\n'
    else:
        header = '# Depender jobs listed first:\n'

    return header + dep_graph_to_text(graph)
