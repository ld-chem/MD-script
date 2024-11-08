from poetry import molgen

# "PAAM.reimage.xml" was got from unwrapped reaction model .xml#

mol1 = molgen.Object("PAAM.reimage.xml",118477, molgen.Shape.none)
gen = molgen.Generators(77.78,77.78,77.78)
gen.addMolecule(mol1,1)

mol2 = molgen.Molecule(1)
mol2.setParticleTypes('S')
gen.addMolecule(mol2,281489)
gen.setMinimumDistance(0.85)
gen.outPutXML("new_box")