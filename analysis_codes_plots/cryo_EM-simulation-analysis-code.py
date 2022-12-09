#	Code to analyse the data from PHDomains fitted inside the cryo-EM map of dynamin assembly
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 7 December 2022


# import the required modules
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.patches as patches
import glob
import seaborn as sns

# Set the fonts to be used in the plots [comment this line if this font is not supported on your machine]
plt.rcParams['font.family'] = 'comfortaa'


# Residue ids of PHDomain
req_resids = range(518,631)


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
for j in req_resids:
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


# Begin plotting of MD simulation data
fig,ax = plt.subplots()
ax.xaxis.label.set_size(12)
ax.yaxis.label.set_size(16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim([req_resids[0],req_resids[-1]])
plt.ylim([-7,45])
l1, =ax.plot(req_resids,meanlist,label = "MD Simulation", color = "C2")
ax.fill_between(req_resids, meanlist-sdlist, meanlist+sdlist,color = "C2",alpha = 0.5)
l3=ax.hlines(y=0, xmin=req_resids[0], xmax=req_resids[-1], colors='grey', linestyles='dashed',label = "Membrane ")
rect1 = patches.Rectangle(xy=(531,-7),width = 6, height = 52,color ="#0075b9", alpha = 0.5,label ="VL-1")
rect4 = patches.Rectangle(xy=(575,-7),width = 6, height = 52,color ="#cc317f", alpha = 0.5,label ="VL-4")
r1 =ax.add_patch(rect1)
r4 =ax.add_patch(rect4)
plt.xlabel("Residue no:")
plt.ylabel("Distance from membrane ("+"$\AA$"+")", fontsize=12,color ="C2")
plt.gca().add_artist(plt.legend(handles = [rect1, rect4],loc="upper right",frameon=False,fontsize=12))
plt.gca().add_artist(plt.legend(handles = [l3],loc="upper left",fontsize=12,frameon=False))
plt.legend(handles=[l1],loc="lower right",fontsize=12,bbox_to_anchor=[1.0, -0.03],frameon=False)

# Add data from cryo-EM analysis in the same plot
ax2  = ax.twinx()

# Get the list of fitted .pdb files in the current directory
fitted_pdbs = glob.glob("./pdbs/*.pdb")

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
	
	# Read the selected .pdb file
	temp_i = open(i,"r").readlines()
	res =0
	# Loop over individual line of .pdb file
	for k in temp_i:
		# Select only C-alpha atom of the residue
		if(k[0:4]=="ATOM" and " CA " in k):
			x, y, z, resi = float(k[30:38]), float(k[38:46]), float(k[46:54]), int(k[22:26])
			
			# Store the distance data of only PHDomain of the entire dynamin molecule
			if(resi>=518 and resi<=630 ):
				# Evaluate distance of the C-alpha atom from the axis of the assembly
				dist_data[res, ctr] = ((y-origin[1])**2 + (z-origin[2])**2)**0.5
				res = res+1
	ctr = ctr+1


# Empty lists for collecting distances of variable loops 1 and 4 from the axis 
# of the assembly from all the fitted PHDomains
v1_collect = []
v4_collect = []

# Matrix to save residue information sorted by relative distance from the axis of the assmebly
ordered_matrix = np.zeros([len(fitted_pdbs),113])

# loop over all the PHDomain data calculated from the individual PHDomain 
for i in range(len(fitted_pdbs)):
	# shifting all the distances in individual PHDomains to observe the orientations on a single plot
	dist_data[:,i] = dist_data[:,i] - np.min(dist_data[:,i])
	
	# plot relative distances of individual residues from the axis of the assembly
	ax2.plot(range(518, 631),dist_data[:,i], color="C7", alpha=0.3)
	ax2.scatter(579,dist_data[(579-518),i], color="#cc317f",marker =".",linewidth=1,s=16)
	ax2.scatter(533,dist_data[(533-518),i], color="#0076ba",marker =".",linewidth=1,s=16)
	
	# Sort the residues by distance from the axis of the assembly
	ordered_matrix[i,:]=np.array(range(518, 631))[np.argsort(dist_data[:,i])]
	
	# Collect the distance of variable loops 1 and 4 from the axis of the assembly
	v1_collect.append(np.min(dist_data[(531-518):(537-518),i]))
	v4_collect.append(np.min(dist_data[(576-518):(583-518),i]))

# Set the plotting parameters
ax2.set_ylabel("Relative Distance from axis of assembly ("+"$\AA$"+")", fontsize=12,color="C4")
plt.tight_layout()
plt.savefig("comparison-cryo_em-simulation_data.png",dpi=500)
plt.close()

# Plot the histogram of distances of variable loops 1 and 4 from the axis 
# of the assembly observed in the fitted PHDomains
counts, bins = np.histogram(v1_collect)
counts = counts/np.sum(counts)
plt.hist(bins[:-1], bins, weights=counts,alpha =0.7,label = "VL1",color="#0076ba" )
counts, bins = np.histogram(v4_collect)
counts = counts/np.sum(counts)
plt.hist(bins[:-1], bins, weights=counts,alpha =0.7,label = "VL4", color = "#cc317f")
plt.legend()
plt.xlabel("Relative Distance from axis of assembly ("+"$\AA$"+")", fontsize=12)
plt.ylabel("Probability", fontsize=12)
plt.tight_layout()
plt.savefig("histogram_VL_1_4.png",dpi=500)
plt.close()



# Initialize a 2D array to store the ranks of individual VLs
loop_matrix = np.zeros([len(fitted_pdbs),4])

# Obtain the ranks of individual VLs using the ordered matrix created above using sorted distances
for i in range(len(fitted_pdbs)):
	temp= []
	for j in range(113):
		if(ordered_matrix[i,j]<=537 and ordered_matrix[i,j]>=531):
			if(1 not in temp):
				temp.append(1)
		if(ordered_matrix[i,j]<=583 and ordered_matrix[i,j]>=576):
			if(4 not in temp):
				temp.append(4)			
		if(ordered_matrix[i,j]<=600 and ordered_matrix[i,j]>=590):
			if(3 not in temp):
				temp.append(3)			
		if(ordered_matrix[i,j]<=560 and ordered_matrix[i,j]>=554):
			if(2 not in temp):
				temp.append(2)
	loop_matrix[i,:] = temp

# Create a map of colors for individual variable loops of PHD
loop_colors = ["#0076ba","#1fb000","#ff9400","#cc317f"]
cmap = mpl.colors.ListedColormap(loop_colors)

# Plot the matrix as an image
plt.matshow(np.transpose(loop_matrix), cmap = cmap,alpha=0.8)
ax = plt.gca()
ax.set_xticks(np.arange(0, 47, 1));
ax.set_yticks(np.arange(0, 5, 1));
ax.set_xticklabels(np.arange(1, 47, 1), fontsize=12)
ax.set_yticklabels(np.arange(1, 5, 1), fontsize=12)
plt.xlabel("PHD id.", fontsize=14)
plt.ylabel("Rank", fontsize=14)
plt.vlines(x=np.arange(0, 46)+0.5, ymin=np.full(4, 0)-0.5, ymax=np.full(4, 4)-0.5, color="black")
plt.hlines(y=np.arange(0, 4)+0.5, xmin=np.full(46, 0)-0.5, xmax=np.full(46, 46)-0.5, color="black")
plt.savefig("loop_rank_matrix.png",dpi=500)
	
