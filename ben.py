import re
import os
import sys
import math
import pydoc
import shutil
import pefile
import numpy as np
import scipy.fftpack
import matplotlib.pyplot as plt

from math import floor
from os.path import join
from tinytree import Tree
from scipy import integrate
from scipy.fftpack import fft

#Period for quad integration
period_length = 10000.0
#Limit for clustering
limit = 0.05
pi_range  = 1

class EntropyNode(Tree):
    """Simple description of entropy block"""

    def __init__(self, off, size, entropy, binary_block, children=None):
        Tree.__init__(self, children)
        self.name = str(off) + str(size) + str(entropy)
        self.offset = off
        self.size = size
        self.entropy = entropy
        binary_block_length = len(binary_block)

        if(binary_block_length < 50):
            return

        lb_length = 0
        rb_length = 0

        if(binary_block_length % 2 == 0):
            rb_length = lb_length = int(binary_block_length / 2)
        else:
            lb_length = int(floor(binary_block_length / 2))
            rb_length = lb_length + 1

        lb_entropy = shannon_entropy(binary_block[0:lb_length])
        rb_entropy = shannon_entropy(
            binary_block[lb_length:binary_block_length])

        # print abs(lb_entropy - rb_entropy), lb_entropy, rb_entropy
        if(abs(lb_entropy - rb_entropy) < 0.1 or lb_entropy == 0 or rb_entropy == 0):
            return

        self.addChild(EntropyNode(self.offset, lb_length,
                                  lb_entropy, binary_block[0:lb_length]))
        self.addChild(EntropyNode(
            self.offset + lb_length, rb_length, rb_entropy, binary_block[lb_length:binary_block_length]))

    def getCrone(self, depth, blocks):
        # print self.offset, self.size, depth, len(blocks)
        if(len(self.children) == 0):
            if depth not in blocks.keys():
                blocks[depth] = []
            blocks[depth].append([self.offset, self.size, self.entropy])
        else:
            self.children[0].getCrone(depth + 1, blocks)
            self.children[1].getCrone(depth + 1, blocks)
        return blocks

    def getDepth(self, previous=0):
        ldepth = previous
        rdepth = previous
        if(len(self.children) != 0):
            ldepth = self.children[0].getDepth(previous + 1)
            rdepth = self.children[1].getDepth(previous + 1)
        if(ldepth > rdepth):
            return ldepth
        else:
            return rdepth

    def dump(self, Nodes=None, Links=None):
        if(not Nodes):
            Nodes = []
        if(not Links):
            Links = []
        if(len(self.children) != 0):
            Nodes.append({"Id": self.offset, "Label": self.name})
            Links.append({"Source": self.name, "Target": self.children[
                         0].name, "Label": str(abs(self.entropy - self.children[0].entropy))})
            Links.append({"Source": self.name, "Target": self.children[
                         1].name, "Label": str(abs(self.entropy - self.children[1].entropy))})
            a, b = self.children[0].dump(Nodes, Links)
            # Nodes.extend(a)
            # Links.extend(b)
            a, b = self.children[1].dump(Nodes, Links)
            # Nodes.extend(a)
            # Links.extend(b)
        return Nodes, Links


def get_execution_block(path_to_file):
    """Return all binary data from executable sections"""
    disfile = path_to_file
    pe = pefile.PE(disfile)

    res = []
    for section in pe.sections:
        if(section.Characteristics & 0x20):
            res.extend(section.get_data())
    return res


def shannon_entropy(binary_data):
    """http://code.activestate.com/recipes/577476-shannon-entropy-calculation/
    :param binary_data:
    """
    freqList = []
    binary_data_length = len(binary_data)
    if(binary_data_length == 0):
        return 8.0
    for b in range(256):
        ctr = 0
        for byte in binary_data:
            if ord(byte) == b:
                ctr += 1
        freqList.append(float(ctr) / binary_data_length)
    ent = 0.0
    for freq in freqList:
        if freq > 0:
            ent = ent + freq * math.log(freq, 2)

    ent = -ent
    return ent


def build_tree(binary_data):
    """
    Build binary tree with entropy nodes
    binary_data - executable section of file
    return - crones of binary tree
    """
    tree = Tree()
    tree.addChild(EntropyNode(0, len(binary_data),
                              shannon_entropy(binary_data), binary_data))
    depth = tree.children[0].getDepth(0)
    blocks = tree.children[0].getCrone(0, {})
    #[self.offset, self.size, self.entropy]
    same_blocks = []
    for index in blocks:
        # print index, len(blocks[index])
        # blocks[index] = map(lambda x : x[0], blocks[index])
        blocks[index].sort(key=lambda x: x[0])
        new_block = [blocks[index][0]]

        current_offset = blocks[index][0][0] + blocks[index][0][1]

        for node_info in blocks[index][1:]:
            if(current_offset == node_info[0]):
                new_block.append(node_info)
                current_offset += node_info[1]
            else:
                same_blocks.append(new_block)
                new_block = []
                new_block.append(node_info)
                current_offset = node_info[0] + node_info[1]
        same_blocks.append(new_block)

    res = []
    for index in same_blocks:
        index.sort(key=lambda x: x[0])

        res.append([index[0][0], index[len(index) - 1][0] +
                    index[len(index) - 1][1], sum(map(lambda x: x[2], index)) / len(index)])

    return res


# FFT
def ff(x, filtered_enb):
    """
    Fast furie function
    filtered_enb - info about entropy blocks (size, entropy, shift)
    return - y
    """
    res = 0.0
    for i in range(len(filtered_enb)):
        m = 1
        if i % 2 != 0:
            m = -1
        block_length = abs(filtered_enb[i][1] - filtered_enb[i][0])
        if block_length == 0: continue
        T = block_length / period_length
        res += m * \
            filtered_enb[i][
                2] * np.sin(x * ((2 * np.pi) / T) + block_length / period_length)
    return res


def quad_fourier(filtered_enb):
    """
    Integrate fft function for fixed period
    filtered_enb - info about entropy blocks (size, entropy, shift)
    """
    x = np.linspace(-pi_range, pi_range)
    y = 0.0
    period_length = 10000.0
    if len(filtered_enb) == 1 and filtered_enb[0][2] == 8.0:
        return 0, 0
    for i in range(len(filtered_enb)):
        m = 1
        if i % 2 != 0:
            m = -1
        block_length = abs(filtered_enb[i][1] - filtered_enb[i][0])
        if block_length == 0: continue
        T = block_length / period_length
        y += m * filtered_enb[i][2] * \
            np.sin(x * ((2 * np.pi) / T) + block_length / period_length)

    yf = fft(y)
    return integrate.quad(lambda a: yf[a], -pi_range, pi_range)



def main():
    """
    arg1 - path to files
    arg2 - path to res dir
    """

    ben_files = []
    corrupted_files = []

    for file_name in os.listdir(sys.argv[1]):
        path_to_file = join(sys.argv[1], file_name)
        if not os.path.isfile(path_to_file): continue
        print "Try to ent-split file", file_name
        try:
            ben_files.append([file_name, build_tree(
                get_execution_block(path_to_file))])
        except pefile.PEFormatError:
            print "Can't process file"
            corrupted_files.append(file_name)

    for ben_file in ben_files:
        ben_file.append(quad_fourier(ben_file[1]))

    ben_sorted_groups = []
    
    for ben_file in ben_files:
        was_found = False
        if(ben_file[2][0] == 0.0): continue
        for group in ben_sorted_groups:
                if abs(group[0][0] - ben_file[2][0]) / group[0][0] < limit and \
                   abs(group[0][1] - ben_file[2][1]) / group[0][1] < limit:
                        was_found = True
                        group[1].append(ben_file[0])
                        break
        if not was_found:
            ben_sorted_groups.append([ben_file[2], [ben_file[0]]])

    unsorted_groups = filter(lambda x: len(x[1]) == 1, ben_sorted_groups)
    ben_sorted_groups = filter(lambda x: len(x[1]) > 1, ben_sorted_groups)

    res_dir = sys.argv[2]

    for group in ben_sorted_groups:
        dir_name = join(res_dir, "[{0}]V={1},E={2}".format(len(group[1]), "%.2f"%group[0][0], "%.2f"%group[0][1]))
        os.mkdir(join(dir_name))
        for file_name in group[1]:        
            shutil.copyfile(join(sys.argv[1], file_name), join(dir_name, os.path.basename(file_name)))


if __name__ == '__main__':
    main()

