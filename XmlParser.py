from io import StringIO
from xml.etree import cElementTree

import numpy as np


# import pandas as pd


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
            # self.data[element.tag] = pd.read_csv(StringIO(element.text),
            #                                    delim_whitespace=True,
            #                                    header=None,
            #                                    ).squeeze("columns").values
            self.data[element.tag] = np.genfromtxt(StringIO(element.text), dtype=None, encoding=None)

