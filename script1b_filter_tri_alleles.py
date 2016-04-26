# Personal Genome Script 1B: filter_tri_alleles.py

# All SNPs in the .map file, rsid index -> rsid
rsids = {}
# Only SNPs that will be processed, mapping rsid index -> rsid
rsids_to_track = {}

# Read in the .map file
with open("samples_hapmap.map") as mapfile:
	count = 0
	for each in mapfile:
		# Make a map of rsid index -> rsid
		rsid = each.strip().split()[1]
		rsids[count] = rsid
		if int(each.strip().split()[3]) >= 0:
			rsids_to_track[count] = rsid
		count = count + 1

num_rsids = len(rsids)
print(len(rsids_to_track))
print(num_rsids)

rsid_allele = {}
# Open the combined .ped file for checking for tri-alleles
with open ("hapmapandme.ped") as pedfile:
	for each in pedfile: 
		pedline = each.strip().split()[6:]
		# For each SNP in the individual's genome
		for i in range(num_rsids):
			# If the SNP is one of the common SNPs we are processing
			if i in rsids_to_track:
				current_snp = rsids_to_track[i]
				# If the SNP has not yet been added, add it to the rsid_allele map
				if current_snp not in rsid_allele:
					rsid_allele[current_snp] = [] 
				
				# If either allele is not in the rsid_allele mapping, add them
				# Don't add the allele if it is unknown (0)
				if pedline[i*2] not in rsid_allele[current_snp] and pedline[i*2] != '0':
					rsid_allele[current_snp].append(pedline[i*2])
				if pedline[i*2+1] not in rsid_allele[current_snp] and pedline[i*2+1] != '0':
					rsid_allele[current_snp].append(pedline[i*2+1])


# Write new .map file and eliminate tri-allelic SNPs 
print(len(rsid_allele))
print("writing to file")
with open("samples_hapmap.map") as infile, open("hapmapandme.map", 'w') as outfile:
	onsnps = 0
	offsnps = 0
	# For each SNP in the mapfile
	for each in infile:
		inline = each.strip().split()
		rsid = inline[1]
		value = int(inline[3])
		# If the rsid is in the rsid_allele map, has more than 2 alleles, and is still a valid SNP
		if rsid in rsid_allele and len(rsid_allele[rsid]) > 2 and value > 0:
			# Negate the value and write the new value to file 
			value = value * -1
			outfile.write(inline[0] + "\t" + inline[1] + "\t" + inline[2] + "\t" + str(value) + "\n")
		else:
			# Write the original line to file
			outfile.write(each)

		# Count the valid or invalid SNPs
		if value > 0:
			onsnps = onsnps + 1
		else:
			offsnps = offsnps + 1
	print(onsnps)
	print(offsnps)
