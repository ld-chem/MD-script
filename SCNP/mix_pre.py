from poetry import molgen

mol2 = molgen.Molecule(1)
mol2.setParticleTypes("S")
mol3 = molgen.Object("out.xml", 100, molgen.Shape.none)
gen=molgen.Generators(10,10,10)
gen.addMolecule(mol3,1)
gen.addMolecule(mol2,750)
gen.setMinimumDistance(0.85)
gen.outPutXML("ABP_chain")
