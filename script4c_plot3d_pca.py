# Personal Genome Script 4C: plot3d_pca.py

import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import pylab
from  matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)

xvals = []
yvals = []
zvals = []
color = []


with open("plink.eigenvec") as infile:
	for each in infile:
		line = each.strip().split()
		xvals.append([float(line[2])])
		yvals.append([float(line[3])])
		zvals.append([float(line[4])])
		if line[0] == "A":
			color.append('green')
		elif line[0] == "B":
			color.append('blue')
		elif line[0] == "C":
			color.append('y')			
		elif not line[1].startswith("NA"):
			color.append('red')
		else:
			color.append('b')
	
	for i in range(len(xvals)):
		ax.scatter(xvals[i], yvals[i], zvals[i], c=color[i])

	plt.show()						
