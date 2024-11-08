# MD-script
My scripts for molecular dynamic simulation of polymer/network.

XmlBuilder.py && __pycache__ && chain_stat.py -> get the index of beads on each chain

XmlParser.py && __pycache__ && density.py -> get the beads density around the center of mass

cm_remake.py -> move the polymer particle to the box center  

catch_bond.py -> catch the bond number of each frame (*.xml)

catch_init.py -> catch the initial beads and change the type of them

dodelete.py && hoomd -> delete the init and cris, rewrite the xml

xml_parser.py && chain.py && cg_topology.py && get_Cross.py -> get the simple closed loop of particles, and determine whether or not chains have crossed the plane

zhuanhuan2.py -> add initial particles, rewrite the xml
