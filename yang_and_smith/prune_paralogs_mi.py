"""
Input: homolog trees
Output: individual orthologs trees

if a tip is longer than the LONG_TIP_CUTOFF
and also long than 10 times its sister, cut it off
This is to fix the leftover trees that frequently has some long tips in it

If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1to1_ORTHOLOGS to False
"""
import newick3, phylo3, os, sys
import importlib
import trim_tree_tips
import tree_utils

OUTPUT_1to1_ORTHOLOGS = True


def get_clusterID(filename):
    return filename.split(".")[0]


def get_front_score(node):
    """
    CJJ: what the bloody hell does this do then
    """
    # print(f'Node is: {node}')
    front_labels = tree_utils.get_front_labels(node)
    # print(f'front_labels:{front_labels}')
    num_labels = len(front_labels)
    # print(f'num_labels: {num_labels}')
    num_taxa = len(set([tree_utils.get_name(i) for i in front_labels]))  # CJJ: tree_utils.get_name(i) just return
    # CJJ: split('.')[0]. So, num_taxa is the number of unique leaf names.
    # print(f'num_taxa: {num_taxa}')
    if num_taxa == num_labels:  # CJJ: i.e. if there are no paralogs
        # print(f'CJJ num_taxa == num_labels', num_taxa)
        return num_taxa
    return -1  # CJJ:  if there are duplicate leaf names (i.e. there are paralogs), return this value


def get_back_score(node, root):
    back_labels = tree_utils.get_back_labels(node, root)
    num_labels = len(back_labels)
    num_taxa = len(set([tree_utils.get_name(i) for i in back_labels]))
    if num_taxa == num_labels:
        return num_taxa
    return -1


def prune(score_tuple, node, root, pp_trees):
    """
    CJJ: and this one?

    score_hashes[node] = (front_score, back_score)
    curroot, done = prune(score_hashes[highest_node], highest_node, curroot, pp_trees)
    """
    if score_tuple[0] > score_tuple[1]:  # prune front  # CJJ: i.e. if more non-repeated taxa in front group
        print("prune front")
        # print(newick3.tostring(node) + ";\n")
        pp_trees.append(node)
        par = node.prune()  # CJJ: par is parent, presumably...
        # CJJ: from phylo3.py:
        # def prune(self):
        #     p = self.parent
        #     if p:
        #         p.remove_child(self)
        #     return p

        if par != None and len(root.leaves()) >= 3:
            par, root = tree_utils.remove_kink(par, root)
        return root, node == root
    else:
        if node != root:  # prune back  # CJJ: i.e. if more non-repeated taxa in back group and it's not the root
            par = node.parent  # par--node<
            par.remove_child(node)
            if par.parent != None:
                par, root = tree_utils.remove_kink(par, root)
        node.prune()
        print("prune back")
        # print(newick3.tostring(root) + ";\n")
        pp_trees.append(root)
        if len(node.leaves()) >= 3:
            node, newroot = tree_utils.remove_kink(node, node)
        else:
            newroot = node
        return newroot, False  # original root was cutoff, not done yet


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("python prune_paralogs_MI.py homoTreeDIR tree_file_ending relative_tip_cutoff absolute_tip_cutoff MIN_TAXA outDIR")
        print("LONG_TIP_CUTOFF is typically same value of the previous LONG_TIP_CUTOFF")
        sys.exit(0)

    inDIR = sys.argv[1] + "/"
    tree_file_ending = sys.argv[2]
    relative_tip_cutoff, absolute_tip_cutoff = float(sys.argv[3]), float(sys.argv[4])
    MIN_TAXA = int(sys.argv[5])
    outDIR = sys.argv[6] + "/"

    for i in os.listdir(inDIR):
        if not i.endswith(tree_file_ending): continue
        print(i)
        with open(inDIR + i, "r") as infile:  # only 1  tree in each file
            intree = newick3.parse(infile.readline())
        curroot = intree
        # print(curroot)
        pp_trees = []

        if get_front_score(curroot) >= MIN_TAXA:  # No need to prune  # CJJ: i.e. no paralogs are more than min taxa
            print("No pruning needed")
            if OUTPUT_1to1_ORTHOLOGS:
                os.system("cp " + inDIR + i + " " + outDIR + get_clusterID(i) + ".1to1ortho.tre")
        else:  # scoring the tree
            going = True
            pp_trees = []

            while going:  # python version of do..while loop
                # print(newick3.tostring(curroot) + ";\n")
                highest = 0
                highest_node = None
                score_hashes = {}  # key is node, value is a tuple (front_score,back_score)
                for node in curroot.iternodes():
                    # front_labels = tree_utils.get_front_labels(node)
                    # num_taxa = len(set([tree_utils.get_name(i) for i in front_labels]))
                    # num_taxa = len(set(front_labels))
                    # print(f'Node {node}, numtaxa: {num_taxa}')

                    front_score = get_front_score(node)  # CJJ: will be > 0 if no paralogs in current clade,
                    # CJJ: or -1 if paralogs
                    back_score = get_back_score(node, curroot)
                    score_hashes[node] = (front_score, back_score)  # CJJ: create dictionary of front and back score
                    # CJJ: for each node
                    if front_score > highest or back_score > highest:
                        highest_node = node  # CJJ: get the node with the greatest number of non-duplicated leaf names (
                        # CJJ: i.e. no paralogs)
                        highest = max(front_score, back_score)  # CJJ: get the number of taxa in the clade with the
                        # CJJ: greatest number of non-duplicated leaf names (i.e. no paralogs)
                if highest >= MIN_TAXA:  # prune
                    # print(f'CJJ curroot is {curroot}')
                    curroot, done = prune(score_hashes[highest_node], highest_node, curroot, pp_trees)
                    # print(f'CJJ curroot is {curroot}')
                    # print(f'CJJ pp_trees is {pp_trees}')
                    # for tree in pp_trees:
                    #     print(newick3.tostring(tree) + ";\n")
                    if done or len(curroot.leaves()) < MIN_TAXA:
                        going = False
                        break
                else:
                    going = False
                    break

        if len(pp_trees) > 0:
            count = 1
            for tree in pp_trees:
                if tree.nchildren == 2:
                    node, tree = trim_tips.remove_kink(tree, tree)
                tree = trim_tips.trim(tree, relative_tip_cutoff, absolute_tip_cutoff)
                # print(f'Tree post-trimming is: {newick3.tostring(tree)};\n')
                if tree != None and len(tree.leaves()) >= MIN_TAXA:
                    # print(f'length: {len(tree.leaves())}\n')
                    # with open(outDIR + get_clusterID(i) + "_MIortho" + str(count) + ".tre", "w") as outfile:
                    with open(outDIR + get_clusterID(i) + ".MIortho" + str(count) + ".tre", "w") as outfile:  # CJJ
                        outfile.write(newick3.tostring(tree) + ";\n")
                    count += 1
