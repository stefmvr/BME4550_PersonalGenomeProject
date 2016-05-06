# Personal Genome Script 4B: plot2d_pca.py

import matplotlib
import matplotlib.pyplot as plt

with open("plink.eigenvec") as infile:
	pcalines = infile.readlines()

	xvals = []
	yvals = []
	code = []

	for each in pcalines:
		pcaline = each.strip().split()
		xvals.append(pcaline[2])
		yvals.append(pcaline[3])

		if pcaline[0] == "A":
			code.append('go')
		elif pcaline[0] == "B":
			code.append('bo')
		elif pcaline[0] == "C":
			code.append('yo')
		else:
			code.append('ro')

	for i in range(len(xvals)):
		plt.plot(xvals[i], yvals[i], code[i])
	
	plt.savefig("pca.png")	
