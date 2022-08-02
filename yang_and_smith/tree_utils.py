#!/usr/bin/env python

# Author: # Author: Yang and Smith
# Modified by: Chris Jackson chris.jackson@rbg.vic.gov.au

from collections import defaultdict

import phylo3
import newick3
import sys


def get_cluster_id(filename):
    """
    Returns first field of tree file name, after splitting by dot/period.

    :param str filename: file name of input tree file
    :return str: first field of tree file name, after splitting by dot/period.
    """

    return filename.split(".")[0]


def get_name(label):
    """
    Given a tip label, return taxon name identifier (first field after splitting by dot/period).

    :param str label: label of a tree tip
    :return:
    """

    return label.split(".")[0]


def get_clusterID(filename):
    """

    :param filename:
    :return:
    """
    return filename.split(".")[0]


def get_front_labels(node):
    """
    Given a node, return a list of front tip labels.

    :param phylo3.Node node: tree object parsed by newick3.parse
    :return list:
    """

    leaves = node.leaves()
    return [i.label for i in leaves]


def get_back_labels(node, root):
    """
    Return taxon names for all child tips OTHER THAN the child tips of the given node

    :param phylo3.Node node: tree object parsed by newick3.parse
    :param phylo3.Node root: tree object parsed by newick3.parse
    :return set:
    """

    all_labels = get_front_labels(root)
    front_labels = get_front_labels(node)
    return set(all_labels) - set(front_labels)  # labels do not repeat


def get_front_names(node):
    """
    Return taxon names for all child tips of the given node. The list may contain identical taxon names.

    :param phylo3.Node node: tree object parsed by newick3.parse
    :return list:
    """

    labels = get_front_labels(node)
    return [get_name(i) for i in labels]


def get_back_names(node, root):

    """
    Given a node, return a list of back tip taxon names. The list may contain identical taxon names.

    :param phylo3.Node node: tree object parsed by newick3.parse
    :param phylo3.Node root: tree object parsed by newick3.parse
    :return list:
    """

    back_labels = get_back_labels(node, root)
    return [get_name(i) for i in back_labels]


def get_front_outgroup_names(node, outgroups):
    """
    Recovers taxon names in tree, and returns a list of the names that are also present in the outgroups list.

    :param phylo3.Node node: tree object parsed by newick3.parse
    :param list outgroups: list of outgroup names recovered from in_and_outgroup_list file
    :return list: a list of taxon names in the provided tree, if they are present in the outgroups list
    """

    names = get_front_names(node)
    return [i for i in names if i in outgroups]


def get_front_ingroup_names(node, ingroups):
    """
    Recovers taxon names in tree, and returns a list of the names that are also present in the ingroups list.

    :param phylo3.Node node: tree object parsed by newick3.parse
    :param list ingroups: list of ingroup names recovered from in_and_outgroup_list file
    :return list: a list of taxon names in the provided tree, if they are present in the ingroups list
    """

    names = get_front_names(node)
    return [i for i in names if i in ingroups]


def remove_kink(node, curroot):
    """
    Smooth the kink created by pruning to prevent creating orphaned tips after pruning twice at the same node.

    :param phylo3.Node node: tree object parsed by newick3.parse
    :param phylo3.Node curroot: tree object parsed by newick3.parse
    :return:
    """

    if node == curroot and curroot.nchildren == 2:
        # move the root away to an adjacent none-tip
        if curroot.children[0].istip:  # the other child is not tip
            curroot = phylo3.reroot(curroot, curroot.children[1])
        else:
            curroot = phylo3.reroot(curroot, curroot.children[0])
    # ---node---< all nodes should have one child only now
    length = node.length + (node.children[0]).length
    par = node.parent
    kink = node
    node = node.children[0]
    # parent--kink---node<
    par.remove_child(kink)
    par.add_child(node)
    node.length = length
    return node, curroot


def pass_boot_filter(node, min_ave_boot):
    """check whether the average bootstrap value pass a cutoff"""
    total = 0.0
    count = 0
    for i in node.iternodes():
        if not i.istip and i.parent != None:
            total += float(i.label)
            count += 1
    if count == 0:  # extracted clades with only two tips
        return True
    boot_average = total / float(count)
    print(boot_average)
    return boot_average >= float(min_ave_boot)


def get_ortho_from_rooted_inclade(inclade,
                                  logger=None):
    """
    Input a rooted tree.
    Cut apart bifurcating nodes when duplicated taxonIDs are detected.

    :param phylo3.Node inclade: tree object parsed by newick3.parse
    :param logging.Logger logger: a logger object
    :return list orthologs: list of orthologs without repeating taxon names
    """

    try:
        assert inclade.nchildren == 2
    except AssertionError:
        logger.error(f'{"[ERROR]:":10} Input clade not properly rooted for clade {newick3.tostring(inclade)}')
        raise

    orthologs = []  # store ortho clades
    clades = [inclade]

    while True:
        newclades = []  # keep track of subclades generated in this round
        for clade in clades:
            num_taxa = len(set(get_front_names(clade)))
            num_tips = len((get_front_labels(clade)))

            if num_taxa == num_tips:  # no taxon repeats
                orthologs.append(clade)
            else:  # has duplicated taxa
                for node in clade.iternodes(order=0):  # PREORDER, root to tip
                    if node.istip:
                        continue
                    # traverse the tree from root to tip
                    try:
                        child0, child1 = node.children[0], node.children[1]
                    except:
                        # print(f'node is {newick3.tostring(node)}')
                        logger.error(f'{"[ERROR]:":10} node.children is: {node.children}')
                        raise

                    name_set0 = set(get_front_names(child0))
                    name_set1 = set(get_front_names(child1))

                    if len(name_set0.intersection(name_set1)) > 0:
                        if node == clade:
                            newclades += [child0, child1]  # break by bifid at the base
                        elif len(name_set0) > len(name_set1):  # cut the side with less taxa
                            node.remove_child(child1)
                            child1.prune()
                            node, clade = remove_kink(node, clade)  # no rerooting here
                            newclades += [clade, child1]
                        else:
                            node.remove_child(child0)
                            child0.prune()
                            node, clade = remove_kink(node, clade)  # no rerooting here
                            newclades += [clade, child0]
                        break

        if not newclades:
            break
        clades = newclades

    return orthologs


def extract_rooted_ingroup_clades(root,
                                  treefile_basename,
                                  ingroups,
                                  outgroups,
                                  min_ingroup_taxa,
                                  logger=None):
    """
    Input a tree with ingroups and at least 1 outgroups.
    Output a list of rooted ingroup clades.

    :param phylo3.Node root: tree object parsed by newick3.parse
    :param str treefile_basename: name of the tree
    :param list ingroups: list of ingroup taxon names
    :param list outgroups: list of outgroup taxon names
    :param int min_ingroup_taxa: minimum number of ingroup taxa required to return a rooted ingroup clade
    :param logging.Logger logger: a logger object
    :return list, collections.defaultdict inclades, inclades_with_fewer_than_min_ingroup_taxa: a list of phylo3.Node
    inclade objects, a dictionary of treename:[inclades with fewer than min_ingroup_taxa]
    """

    inclades_list = []
    inclades_with_fewer_than_min_ingroup_taxa_list = []

    while True:
        max_score, direction, max_node = 0, "", None
        for node in root.iternodes():
            front, back = 0, 0
            front_names_set = set(get_front_names(node))

            for name in front_names_set:
                if name in outgroups:
                    front = -1
                    break
                elif name in ingroups:
                    front += 1
                else:
                    sys.exit("Check taxonID " + name)
            back_names_set = set(get_back_names(node, root))

            for name in back_names_set:
                if name in outgroups:
                    back = -1
                    break
                elif name in ingroups:
                    back += 1
                else:
                    sys.exit("Check taxonID " + name)

            if front > max_score:
                max_score, direction, max_node = front, "front", node
            if back > max_score:
                max_score, direction, max_node = back, "back", node

        if max_score >= min_ingroup_taxa:
            if direction == "front":
                inclades_list.append(max_node)
                kink = max_node.prune()
                if len(root.leaves()) > 3:
                    newnode, root = remove_kink(kink, root)
                else:
                    break

            elif direction == "back":
                par = max_node.parent
                par.remove_child(max_node)
                max_node.prune()
                inclades_list.append(phylo3.reroot(root, par))  # flip direction
                if len(max_node.leaves()) > 3:
                    max_node, root = remove_kink(max_node, max_node)
                else:
                    break

        elif max_node:
            logger.debug(f'Clade {newick3.tostring(max_node)} from tree {treefile_basename} contains fewer than the '
                         f'min_ingroup_taxa value of {min_ingroup_taxa}. Skipping clade.')
            inclades_with_fewer_than_min_ingroup_taxa_list.append(newick3.tostring(max_node))
        else:
            break

    return inclades_list, inclades_with_fewer_than_min_ingroup_taxa_list
