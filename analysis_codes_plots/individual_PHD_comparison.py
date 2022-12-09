#	Code to plot the individual PHDomain orientation along with that observed in simulation 
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 7 December 2022

# Create an empty directory named "individual-phd-comparison" before executing this script


# import the required modules
import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as patches
import glob

# Set the fonts to be used in the plots [comment this line if this font is not supported on your machine]
plt.rcParams['font.family'] = 'comfortaa'

# Residue ids of PHDomain
req_resi = range(518,631)


# Read C_alpha distances from the data derived from all-atom MD Simulation
# Each data file correspond to individual residue of PHDomain 
# pos-10.dat contains the distance of 518 (10+508) residue C_alpha atom from 
# the plane of lipid bilayer 

# get trajectory length and time frames from the first file in the simulation-residues-z-height folder
traj = open("./simulation-residues-z-height/pos-10.dat").readlines()
data = np.empty(shape = [len(traj),114])
for k in range(len(traj)):
	data[k,0] =  float(traj[k].split()[0])
	
# Loop over all the residues of PHDomain
counter  = 1
for j in req_resi:
	traj = open("./simulation-residues-z-height/pos-"+str(j-508)+".dat").readlines()
	for k in range(len(traj)):
		data[k,counter] =float(traj[k].split()[1])
	counter  = counter +1

# Evaluate Mean and standard deviation of positions observed over the last 
# few nanoseconds of the MD simulation
meanlist = []
sdlist = []
for k in range(1,114):
	meanlist = meanlist + [np.mean(data[-20:,k])]
	sdlist = sdlist + [np.std(data[-20:,k])]
meanlist, sdlist = np.array(meanlist), np.array(sdlist)


# Get the list of fitted .pdb files in the current directory
fitted_pdbs = glob.glob("pdbs/*.pdb")

# Initialize an empty 2D array to store the distance data of individual residue(C_alpha atom) 
# of individual PHDs from the axis of the dynamin assembly on the membrane
dist_data = np.zeros([113,len(fitted_pdbs)])


ctr = 0
# address of the central voxel in the cryo-EM map
central_voxel = np.array([283, 283, 283])

# Convert discrete voxel into real units (angstroms)
origin = central_voxel * 1.07


# Loop over all the .pdb files
for i in fitted_pdbs:
	i_x, i_y, i_z, dist, res_no = [], [], [], [], []
	temp_i = open(i,"r").readlines()
	res =0
	for k in temp_i:
		if(k[0:4]=="ATOM" and " CA " in k):
			x, y, z, resi = float(k[30:38]), float(k[38:46]), float(k[46:54]), int(k[22:26])
			if(resi>=518 and resi<=630 ):
				dist_data[res, ctr] = ((y-origin[1])**2 + (z-origin[2])**2)**0.5
				res = res+1
	ctr = ctr +1

# Begin plotting of cryo-EM analysis data
for i in range(len(fitted_pdbs)):
	fig,ax = plt.subplots()
	ax.xaxis.label.set_size(12)
	ax.yaxis.label.set_size(14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.xlim([req_resi[0],req_resi[-1]])
	plt.ylim([-7,45])
	dist_data[:,i] = dist_data[:,i] - np.min(dist_data[:,i])
	l1,=ax.plot(range(518, 631),dist_data[:,i], alpha=0.8,color ="C4",label = "Cryo EM data")
	l3=ax.hlines(y=0, xmin=req_resi[0], xmax=req_resi[-1], colors='grey', linestyles='dashed',label = "Membrane ")
	rect1 = patches.Rectangle(xy=(531,-7),width = 6, height = 52,color ="#0075b9", alpha = 0.5,label ="VL-1")
	rect4 = patches.Rectangle(xy=(575,-7),width = 6, height = 52,color ="#cc317f", alpha = 0.5,label ="VL-4")
	r1 =ax.add_patch(rect1)
	r3 =ax.add_patch(rect4)
	plt.xlabel("Residue no:")
	plt.ylabel("Relative Distance from axis of assembly ("+"$\AA$"+")", fontsize=12,color = "C4")	
	plt.gca().add_artist(plt.legend(handles = [rect1, rect4],loc="upper right",frameon=False,fontsize=12))
	plt.gca().add_artist(plt.legend(handles = [l3],loc="upper left",fontsize=12,frameon=False))

	#Add data from MD simulation data in the same plot
	ax2  = ax.twinx()
	l2,=ax2.plot(req_resi,meanlist,label = "MD Simulation", color = "C2")
	ax2.fill_between(req_resi, meanlist-sdlist, meanlist+sdlist,color ="C2",alpha = 0.5)		
	ax2.set_ylabel("Distance from membrane ("+"$\AA$"+")", fontsize=12,color = "C2")
	plt.legend(handles=[l1,l2],loc="lower right",fontsize=12,bbox_to_anchor=[1.0, -0.03],frameon=False)
	plt.tight_layout()
	plt.savefig("./individual-phd-comparison/"+str(i)+".png",dpi=300)
	plt.close()

