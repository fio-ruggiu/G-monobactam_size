#!/usr/bin/env python

import argparse
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import numpy as np
import pandas as pd
import itertools
import math
import re
#from matplotlib import pyplot
#from mpl_toolkits.mplot3d import Axes3D

# Set up argument parser
def read_args():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--input_file', required=True, type=str, help='MOL2 with conformers on which to calculate min dimensions')
	parser.add_argument('--output_file', default='min_dimensions.txt', required=False, type=str, help='File name for outputs')
	parser.add_argument('--min_dim1', default=10, required=False, type=float, help='First minimal dimension')
	parser.add_argument('--min_dim2', default=14, required=False, type=float, help='Second minimal dimension')
	parser.add_argument('--rot_angle', default=5, required=False, type=float, help='Angle increment for rotation of object')
	parser.add_argument('--use_convex_hull', default=False, required=False, type=bool, help='Boolean to whether to use convex hull or not. If set to False, then all coordinates will be used. Defaults to False')
	
	args = parser.parse_args()
	return args

def vanderWaals_radii(elt_symbol): #van der Waals radii in A from Bondi,J.Phys.Chem.,68,441, 1964. B vdW radii is the estimate according to Pauling. Other elements are assigned van der Waals radii of 2.00A.
	vdW_radii = {
			"H"	: 1.20,
			"C"	: 1.70,
			"N"	: 1.55,
			"O"	: 1.52,
			"F"	: 1.47,
			"Cl"	: 1.75,
			"Br"	: 1.85,
			"I"	: 1.98,
			"S"	: 1.80,
			"P"	: 1.80,
			"B"	: 1.65,
			"Ag"	: 1.72,
			"Ar"	: 1.88,
			"As"	: 1.85,
			"Au"	: 1.66,	
			"Cd"	: 1.58,		
			"Cu"	: 1.40,		
			"Ga"	: 1.87,		
			"He"	: 1.40,
			"Hg"	: 1.55,		
			"In"	: 1.93,
			"K"	: 2.75,
			"Kr"	: 2.02,
			"Li"	: 1.82,
			"Mg"	: 1.73,
			"Na"	: 2.27,
			"Ne"	: 1.54,
			"Ni"	: 1.63,
			"Pb"	: 2.02,
			"Pd"	: 1.63,
			"Pt"	: 1.72,
			"Se"	: 1.90,
			"Si"	: 2.10,
			"Sn"	: 2.17,
			"Te"	: 2.06,
			"Tl"	: 1.96,
			"U"	: 1.86,
			"Xe"	: 2.16,
			"Zn"	: 1.39,
		}
	
	if (vdW_radii[elt_symbol]): return vdW_radii[elt_symbol]
	else: return 2.00


def RetrieveMol2Block(fileLikeObject, delimiter="@<TRIPOS>MOLECULE"):
    """generator which retrieves one mol2 block at a time
    """
    mol2 = []
    for line in fileLikeObject:
        if line.startswith(delimiter) and mol2:
            yield "".join(mol2)
            mol2 = []
        mol2.append(line)
    if mol2:
        yield "".join(mol2)

def RetrieveMol2Properties(mol, mol_block): 
    """Retrieves mol2 properties and input them in rdkit mol object
    """
    for result in re.findall('@<TRIPOS>PROPERTY_DATA(.*?)#', mol_block, re.S):
    	for line in result.split('\n'):
    		prop=re.sub('\s+','',line).split('|')
    		if (len(prop)==2):
    			mol.SetProp(prop[0],prop[1])
    			#print(mol.GetProp(prop[0]))

def rotate_x(convex_hull,angle): #rotate around y axis by angle (rad)
   new_convex_hull=np.empty((0,3), float)
   for point in convex_hull:
   	new_convex_hull=np.append(new_convex_hull, values=[[point[0] , point[1]*math.cos(angle)-point[2]*math.sin(angle) , point[1]*math.sin(angle)+point[2]*math.cos(angle) ]],axis=0) 	
    	#new_convex_hull=np.append(new_convex_hull, values=np.dot(np.array([[0, math.cos(math.radians(rot_angle_x)), -math.sin(math.radians(rot_angle_x))],[0, math.sin(math.radians(rot_angle_x)), -math.cos(math.radians(rot_angle_x))]]),point))
   return new_convex_hull

def rotate_y(convex_hull,angle): #rotate around y axis
   new_convex_hull=np.empty((0,3), float)
   for point in convex_hull:
    	new_convex_hull=np.append(new_convex_hull, values=[[ point[0]*math.cos(angle)-point[2]*math.sin(angle) , point[1] , point[0]*math.sin(angle)+point[2]*math.cos(angle)]],axis=0)
    	#new_convex_hull=np.append(new_convex_hull, values=np.dot(np.array([[0, math.cos(math.radians(rot_angle_x)), -math.sin(math.radians(rot_angle_x))],[0, math.sin(math.radians(rot_angle_x)), -math.cos(math.radians(rot_angle_x))]]),point))
   return new_convex_hull

def rotate_z(convex_hull,angle): #rotate around z axis
   new_convex_hull=np.empty((0,3), float)
   for point in convex_hull:
    	new_convex_hull=np.append(new_convex_hull, values=[[ point[0]*math.cos(angle)-point[1]*math.sin(angle) , point[0]*math.sin(angle)+point[1]*math.cos(angle) , point[2]]],axis=0)
    	new_convex_hull=np.append(new_convex_hull, values=np.dot(np.array([[0, math.cos(math.radians(rot_angle_x)), -math.sin(math.radians(rot_angle_x))],[0, math.sin(math.radians(rot_angle_x)), -math.cos(math.radians(rot_angle_x))]]),point))
   return new_convex_hull

#calculate angle between dipole moment vector and z-axis:
def calc_theta(dip_mom):
	theta=math.degrees(math.acos(np.dot(dip_mom,np.array([0,0,1]))/(np.linalg.norm(dip_mom)*np.linalg.norm(np.array([0,0,1])))))
	return theta
	
#Project convex hull points onto x,y plane and check measurements
def max_dim_calc(convex_hull,min_dim1,min_dim2,vdW_radii):
	min_dim_bool=False
	x_max=0 #Maximum distance between two atoms in the x projection
	y_max=0 #Maximum distance between two atoms in the y projection
	z_max=0 #Maximum distance between two atoms in the z projection

	indices_comb=list((i,j) for ((i,_),(j,_)) in itertools.combinations(enumerate(convex_hull), 2))
	#for atm_comb in list(itertools.combinations(convex_hull,2)):
	for i,j in indices_comb:
		#print(i,j,convex_hull[i][0],convex_hull[j][0],vdW_radii[i],vdW_radii[j])
		diff_x=abs(convex_hull[i][0]-convex_hull[j][0]) +vdW_radii[i]+vdW_radii[j]
		diff_y=abs(convex_hull[i][1]-convex_hull[j][1]) +vdW_radii[i]+vdW_radii[j]
		diff_z=abs(convex_hull[i][2]-convex_hull[j][2]) +vdW_radii[i]+vdW_radii[j]
		if(x_max<diff_x):x_max=diff_x
		if(y_max<diff_y):y_max=diff_y
		if(z_max<diff_z):z_max=diff_z
	if((x_max<=min_dim1)and(y_max<=min_dim2)) or ((x_max<=min_dim2)and(y_max<=min_dim1)):
		min_dim_bool=True
	
	return min_dim_bool,x_max,y_max,z_max
	
if __name__ == '__main__':
	args = read_args()

	#Read in molecular file and loop over each molecule
	#confSupp=Chem.SDMolSupplier(args.input_file, removeHs=False)
	#for mol in confSupp:
	outfile = open(args.output_file, 'w')
	outfile.write("Conf_nb,dipole_moment_magnitude,Bool_min_dim,min_area,min_theta,min_area_theta, min_area_x,min_area_y,min_area_z,min_theta_area,min_theta_x,min_theta_y,min_theta_z\n")
	outfile.close()
	conf_nb=0
	with open(args.input_file) as fi:
		for mol2 in RetrieveMol2Block(fi):
			mol = Chem.MolFromMol2Block(mol2, removeHs=False)
			RetrieveMol2Properties(mol, mol2) 
			conf_nb+=1
			print("conformer number: "+str(conf_nb))
			#Center molecule on x,y,z origin
			Chem.rdMolTransforms.CanonicalizeConformer(mol.GetConformer(),ignoreHs=False)

			#Translate x,y,z coordinates of molecule object to numpy array
			coord=np.array(mol.GetConformer().GetPositions())
			#coord=np.array([[0,1,0],[0,-1,0]])
			
			#Get each atoms' partial charges
			partial_charges = np.array([float(a.GetProp("_TriposPartialCharge")) for a in mol.GetAtoms()])
			#partial_charges = np.array([-0.5,0.5])
			
			#Get each atoms' vdW radii
			#vdW_radii_orig = np.array([float(Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum())) for a in mol.GetAtoms()])
			vdW_radii_orig = np.array([float(vanderWaals_radii(a.GetSymbol())) for a in mol.GetAtoms()])
			
			#Calculate dipole moment based on coordinates and partial charges. Note: Center of molecule corresponds to origin
			partial_charges_diag = np.diag(partial_charges)
			dip_mom=np.sum(np.dot(partial_charges_diag,coord), axis=0)
			dip_mom_magnitude=np.linalg.norm(dip_mom)
			#print(partial_charges[0], coord[0], partial_charges_diag[0])
			#theta=math.degrees(math.acos(np.dot(dip_mom,np.array([0,0,1]))/(np.linalg.norm(dip_mom)*np.linalg.norm(np.array([0,0,1])))))
			#print(dip_mom, np.linalg.norm(dip_mom), theta)
			
			# calculate convex hull
			if (args.use_convex_hull==True):
				from scipy.spatial import ConvexHull
				convex_hull_orig=ConvexHull(coord)
				convex_hull=np.array(convex_hull_orig.points[convex_hull_orig.vertices])
				vdW_radii=np.array(vdW_radii_orig[convex_hull_orig.vertices])
				
			else: 
				convex_hull=coord
				vdW_radii=vdW_radii_orig
			
			#Initialize
			min_area_values=["False",99999999,360,999999,999999,999999] #array with dimensions corresponding to minimal rectangle area [min_area,x,y,z]
			min_dim_bool=0
			theta_min=[99999999,360,999999,999999,999999]
			
			rot_angle_x=0
			rot_angle_y=0
			while(rot_angle_x<359):
				#Rotate convex hull around x-axis
				dip_mom_x=rotate_x(np.reshape(dip_mom, (1,3)),math.radians(rot_angle_x))				
				convex_hull_x=rotate_x(convex_hull,math.radians(rot_angle_x))
				
				#Do all combinations of points in the convex hull and check to find max distances in each direction (x,y)
				min_dim,x_max,y_max,z_max=max_dim_calc(convex_hull_x,args.min_dim1,args.min_dim2,vdW_radii)
				theta=calc_theta(dip_mom_x)
				
				if(min_dim==True): 
						min_dim_bool=1
						if(theta_min[1]>theta):theta_min=[x_max*y_max,theta,x_max,y_max,z_max]
				if(min_area_values[1]>x_max*y_max): min_area_values=[str(min_dim),x_max*y_max,theta,x_max,y_max,z_max]
				#print(min_dim_bool,rot_angle_x,rot_angle_y,x_max*y_max,x_max,y_max,z_max,theta)
				#if(min_dim_bool): break
				rot_angle_y=0
				while(rot_angle_y<359):
					#Rotate convex hull around y-axis		
					convex_hull_y=rotate_y(convex_hull_x,math.radians(rot_angle_y))
					dip_mom_y=rotate_y(dip_mom_x,math.radians(rot_angle_y))
					min_dim,x_max,y_max,z_max=max_dim_calc(convex_hull_y,args.min_dim1,args.min_dim2,vdW_radii)
					theta=calc_theta(dip_mom_y)
					
					if(min_dim==True): 
						min_dim_bool=1
						if(theta_min[1]>theta):theta_min=[x_max*y_max,theta,x_max,y_max,z_max]
					if(min_area_values[1]>x_max*y_max): min_area_values=[str(min_dim),x_max*y_max,theta,x_max,y_max,z_max]
					#print(min_dim_bool,rot_angle_x,rot_angle_y,x_max*y_max,x_max,y_max,z_max,theta)
					#if(min_dim_bool): break
					rot_angle_y+=args.rot_angle
				#if(min_dim_bool): break
				rot_angle_x+=args.rot_angle
			
			outfile = open(args.output_file, 'a')
			outfile.write("{:d},{:.2f},{:d},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n".format(conf_nb,dip_mom_magnitude,min_dim_bool,min_area_values[1],theta_min[1],min_area_values[2],min_area_values[3],min_area_values[4],min_area_values[5],theta_min[0],theta_min[2],theta_min[3],theta_min[4]))
			outfile.close()
				

