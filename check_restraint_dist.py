#!/usr/bin/env python

import argparse
import math
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles

# Set up argument parser
def read_args():
	parser = argparse.ArgumentParser(description='Merge conformer with NMR structure, align atoms to reference molecule and extract distances in conformers of each restraint in each conformer')
	parser.add_argument('--conf_file', required=True, type=str, help='SDF with conformers to be compared to reference molecule. Note: Should include hydrogen positions.')
	parser.add_argument('--output_file', default='output.sdf', required=False, type=str, help='Output SDF with restraint distances, Boltzmann population, and flag if both restraint within cut-off for each conformer')
	parser.add_argument('--eflag', default='E', required=False, type=str, help='SDF field name for energies/probabilities of conformers to be extracted. Indicate None if there is no energy field.')
	parser.add_argument('--cutoff', default=7, required=False, type=int, help='Cut-off to be applied to restraints. By default 7A')
	args = parser.parse_args()
	return args

def conf_boltzmann_eval(mol, conformers, energies):
	k_B = 0.0019872041; #in kcal mol-1 K-1 as MOE gives energies in kcal mol-1
	T = 298.15; #Room temperature in K.
	k_BT = k_B*T
	sum_exp_EkBT=0
	sum_Aexp_EkBT=0
	for i in range(0,len(conformers),1):
		sum_exp_EkBT += math.exp(-energies[i]/k_BT);
	
	sum_Aexp_EkBT += imhb_conf*math.exp(-energies[i]/k_BT);	
	Bavg = sum_Aexp_EkBT/sum_exp_EkBT;
	
	return 
	
if __name__ == '__main__':

	args = read_args()
	print(args)
	refMol=Chem.MolFromSmiles('[H]N([H])c1nc(c([H])s1)C(=N\OC([H])(C([O-])=O)C([H])([H])Oc1c([H:5])c([H:2])c(c([H:4])c1[H:1])C(N([H])C1([H])C([H])([H])C([H])([H])[N+]([H])([H])C1([H])[H])=[N+]([H])[H])\C(=O)N([H])C1([H:28])C(=O)N(C1([H])C([H])([H])n1nc([H])nc1[H:38])S([O-])(=O)=O', sanitize=False) #note that molecule should already be rdkit sanitized in this SMILES
	#Create dictionary between mapped atoms and atom numbers
	atom_dict={}
	for atom in refMol.GetAtoms():
		if(atom.GetAtomMapNum()):
			atom_dict[atom.GetAtomMapNum()]=atom.GetIdx()
			#atom.SetAtomMapNum(0)

	#read in conformers
	confSupp=Chem.SDMolSupplier(args.conf_file, removeHs=False)
	conformers=[]
	energies=[]
	conf_nb=[]
	conf_count=0
	for m in confSupp:
		conf_count+=1
		if(args.eflag=='None'):
			energies.append("1")
			conformers.append(m)
			conf_nb.append(conf_count)
		else:
			try:
				energies.append(float(m.GetProp(args.eflag)))
				conformers.append(m)
				conf_nb.append(conf_count)
			except:
				print("Can't read energy for conformer: "+str(conf_count)+" - ignoring conformer!!!")

	#Match reference onto conformers to get atom mapping
	ssm=conformers[0].GetSubstructMatches(refMol)[0]
	
   #Define complete atom map for renumbering
	for i in range(0,conformers[0].GetNumAtoms()):
		if i not in ssm:
			ssm.append(i)
			
	#Calculate Boltzmann averaged
	k_B = 0.0019872041; #in kcal mol-1 K-1 as MOE gives energies in kcal mol-1
	T = 298.15; #Room temperature in K.
	k_BT = k_B*T
	sum_exp_EkBT=0
	for i in range(0,len(conformers),1):
		sum_exp_EkBT += math.exp(-energies[i]/k_BT);
		
	writer = Chem.SDWriter(args.output_file)
	#Loop over the conformer molecules and renumber them
	for i in range(0,len(conformers)):
		m=rdmolops.RenumberAtoms(conformers[i],ssm) #renumber molecules to match atom numbering in RefMol
		
		#calculate distances of restraints and indicate if both average restraints =<5Angstrom in the conformer
		restraint_1_dist=(Chem.rdMolTransforms.GetBondLength(m.GetConformer(),atom_dict[2],atom_dict[38])+Chem.rdMolTransforms.GetBondLength(m.GetConformer(),atom_dict[4],atom_dict[38]))/2
		restraint_2_dist=(Chem.rdMolTransforms.GetBondLength(m.GetConformer(),atom_dict[1],atom_dict[28])+Chem.rdMolTransforms.GetBondLength(m.GetConformer(),atom_dict[5],atom_dict[28]))/2
		
		m.SetProp("restraint_1_dist",str(restraint_1_dist))
		m.SetProp("restraint_2_dist",str(restraint_2_dist))
		if((restraint_1_dist<=args.cutoff)and(restraint_2_dist<=args.cutoff)): m.SetProp("NMR_conf","True")
		else: m.SetProp("NMR_conf","False")
		
		m.SetProp("E",str(energies[i]))
		m.SetProp("glob",conformers[i].GetProp("glob"))
		m.SetProp("rgyr",conformers[i].GetProp("rgyr"))
		Bavg=100*math.exp(-energies[i]/k_BT)/sum_exp_EkBT 
		m.SetProp("Boltzmann_pop(%)",str(Bavg))
		
		writer.write(m)
	writer.close()
