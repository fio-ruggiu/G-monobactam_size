#Analyzing gram-negative entry of compound 1
#Author: Fiorella Ruggiu, Novartis Institutes for BioMedical Research, ruggiu.fiorella@gmail.com
#Last reviewed: 2018/10/05

#Generate conformers for compound1
moebatch -exec "file_confsearch['compound1_pKa.sdf','sdf','compound1_MOE_e80_100_Ele.sdf','sdf',[maxconf:100,solDielectric:80,eleEnable:1]]"
moebatch -exec "file_confsearch['compound1_pKa.sdf','sdf','compound1_MOE_e80_100_noEle.sdf','sdf',[maxconf:100,solDielectric:80,eleEnable:0]]"

#Calculate restraints corresponding to restraint 1 and 2 (RDKit 2018.03.4.0 and Anaconda 3.5-4.0.0)
./check_restraint_dist.py --conf_file compound1_MOE_e80_100_Ele.sdf --output_file compound1_MOE_e80_100_Ele_restraints.sdf
./check_restraint_dist.py --conf_file compound1_MOE_e80_100_noEle.sdf --output_file compound1_MOE_e80_100_noEle_restraints.sdf

#E.coli OmpF + ampicillin crystal structure (PDB ID: 4gcp) and superpose conformers
##Open 4gcp PDB file
##Compute -> Protonate 3D
##Open conformer SDF into mdb database
##Select lactam+amide motif from ampicillin and do Compute->Molecule->Superpose in database menu
##Surface of receptor calculated in main window Surface->Surfaces and Maps... with options:
###Name: 4GCP
###Atoms: Receptor Atoms
###Near: Receptor Atoms
###Color: Electrostatics
##Calculate Amber10:EHT charges (Compute->Molecule->Partial Charges) and export stable and restraint-compliant conformers (NMR_conf==True && Boltzmann_pop(%) >= 2) to mol2 format (compound1_conf1.mol2 and compound1_conf2-4.mol2)
##Note: to draw dipole use Drawdipole[1] in svl command line.
##Note: to calculate dipole use Compute->Descriptors->AM1_dipole/dipole in database menu
##PDB and MOE file with 4 stable conformer available in folder (4gcp_compound1_conformer_superposed.pdb/moe)

#Minimal dimensions and dipole moment angle calculations (RDKit 2018.03.4.0, numpy 1.15.1, pandas 0.23.4, scipy 1.1.0 and Anaconda 3.5-4.0.0)
./min_dim_mol.py --input_file compound1_conf1.mol2 --output_file compound1_conf1_min_dim.txt
./min_dim_mol.py --input_file compound1_conf2-4.mol2 --output_file compound1_conf2-4_min_dim.txt
##Note: min_area and min_area_theta are the parameters reported in the publication
##Note edited compound1_conf2-4_min_dim.txt to have corresponding conf numbers
