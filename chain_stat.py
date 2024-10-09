import networkx as nx
from sys import argv
import numpy as np
import re
from io import StringIO
from xml.etree import cElementTree


### Universal script for chain index statistic. ###
### Network also works. ###


def control_in(control_file):
    pass


class Box(object):
    def __init__(self):
        self.xy = 0
        self.xz = 0
        self.yz = 0
        return

    def update(self, dic):
        self.__dict__.update(dic)


class XmlParser(object):
    def __init__(self, filename, needed=None):
        tree = cElementTree.ElementTree(file=filename)
        root = tree.getroot()
        self.box = Box()
        self.data = {}
        needed = [] if needed is None else needed
        for key in root[0].attrib:
            self.__dict__[key] = int(root[0].attrib[key])
        for element in root[0]:
            if element.tag == 'box':
                self.box.update(element.attrib)
                continue
            if (len(needed) > 0) and (element.tag not in needed):
                continue
            if element.tag == 'reaction':
                self.data['reaction'] = []
                reaction_list = element.text.strip().split('\n')
                while '' in reaction_list:
                    reaction_list.remove('')
                for l in reaction_list:
                    r = re.split(r'\s+', l)
                    while '' in r:
                        r.remove('')
                    r[1:] = [int(_) for _ in r[1:]]
                    self.data['reaction'].append(r)
                continue
            if element.tag == 'template':
                self.data['template'] = eval('{%s}' % element.text)
                continue
            if len(element.text.strip()) > 0:
                self.data[element.tag] = np.genfromtxt(StringIO(element.text), dtype=None, encoding=None)

top = nx.Graph()
xml = XmlParser(argv[1])
position_ = xml.data['position']
types_ = xml.data['type']
#types_ = types_[types_ == 'A']
for i,t in enumerate(xml.data['type']):
    top.add_node(i)
for b in xml.data['bond']:
#    if b[1] == "A" and b[2] == "A":
    top.add_edge(b[1],b[2],bondtype=b[0])

visited = set()

def explore_all_connected_nodes(start_node, f):
    queue = [start_node]
    all_related_nodes = set()
    while queue:
        node = queue.pop()
        if node in visited:
            continue
        visited.add(node)
        all_related_nodes.add(node)
        neighbors = list(top.neighbors(node))

#        f.write(f"Node {node} is connected to: {neighbors}\n")

        for neighbor in neighbors:
            if neighbor not in visited:
                queue.append(neighbor)

    f.write(f"All nodes connected to node {start_node}: {sorted(all_related_nodes)}\n")

with open('output.txt', 'w') as f:
    starting_node = 1
    explore_all_connected_nodes(starting_node, f)

    for node in top.nodes():
        if node not in visited:
            explore_all_connected_nodes(node,f)
