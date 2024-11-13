# MD-script
My scripts for molecular dynamic simulation of polymer/network.

XmlBuilder.py && __pycache__ && chain_stat.py && chain_stat_transfer.py -> 
    get the index of beads on each chain, and change the form for ovito

XmlParser.py && __pycache__ && density.py -> get the beads density around the center of mass

cm_remake.py -> move the polymer particle to the box center  

catch_bond.py -> catch the bond number of each frame (*.xml)

catch_init.py -> catch the initial beads and change the type of them

dodelete.py && hoomd -> delete the init and cris, rewrite the xml

xml_parser.py && chain.py && cg_topology.py && get_Cross.py -> get the simple closed loop of particles, and determine whether or not chains have crossed the plane

zhuanhuan2.py -> add initial particles, rewrite the xml

SCNP -> make model of single chain nanoparticles

PAAM -> simulate the hydrogel reaction in situ

ddmonomer.py -> delete the monomers in system

1a.py && 3a.py && tension.py && tensile.py && draw_tensile.py -> uniaxial and triaxial tensile tests and stress-strain curve drawing

resize.py -> additional solvents were added to the model to simulate the swelling process of the hydrogel, then nvt 

monomer_calculate.py -> calculate the monomer conversion

catch_Rg.py -> calculate the variation of Rg during uniaxial tensile

rexml.py -> delete the beads that do not need anymore, rewrite one .xml
