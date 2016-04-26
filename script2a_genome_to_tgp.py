# Personal Genome Script 2B: genome_to_tgp.py

# .map file dictionary, maps the rsid index -> rsid
map_dict = {}
# .map file dictionary, maps the rsid -> rsid index (limits genome reads to only necessary SNPs)
rsid_map_dict = {}
with open("tgp_plink.map") as mapfile:
	count = 0
	for each in mapfile:
		mapline = each.strip().split()
		map_dict[count] = mapline[1]
		rsid_map_dict[mapline[1]] = count
		count = count + 1
num_genes = len(map_dict)
print(num_genes)

# Mapping of rsid -> the number of sample individuals that include this SNP
gen_count = {}
# Possible chromosome suffixes for each chromosome file of an individual
suffixes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20','21','22','X']
# Helper Function to process the imputed genomes for each of the samples
def process_gen(dirname, gid, fileid):
	gen_dict = {}
	
	# Create a new .ped file for the sample's genome
	with open(gid+"_tgp.ped", 'w') as outfile:
		# For each of the chromosomes
		for suffix in suffixes:
			# Open the file associated with that chromosome
			filename = dirname + "/id_" + fileid + "_chr" + suffix + ".23andme.txt"
			with open(filename) as genfile:
				# For each RSID in the file, map the RSID to the genotype call
				for each in genfile:
					genline = each.strip().split()
					if genline[0] in rsid_map_dict:
						gen_dict[genline[0]] = genline[3]
					
						# Count the occurance of each SNP to identify common SNPs
						if genline[0] not in gen_count:
							gen_count[genline[0]] = gid
						elif gid not in gen_count[genline[0]]:
							gen_count[genline[0]] = gen_count[genline[0]] + gid
		# Write the first 6 columns of the .ped file
		outfile.write(gid + "\t" + gid + "\t0\t0\t0\t-9\t")	

		count_valid = 0
		# Iterate through all of the SNPs in the .map file
		for i in range(num_genes):
			rsid = map_dict[i]
			# If the RSID is present and the call has both alleles
			# (The X chromosome is already homozygous in the imputed data)
			if rsid in gen_dict and len(gen_dict[rsid]) == 2:
				# Write the call to the output file
				call = gen_dict[rsid]
				outfile.write(call[0] + " " + call[1] + "\t")
				count_valid = count_valid + 1
			else:
				outfile.write("0 0\t")
		outfile.write("\n")
		print(count_valid)

# Make calls to process each imputed genome
process_gen("Genome_A_Imputed", "A", "48821H913")
process_gen("Genome_B_Imputed", "B", "4w62L2v54")
process_gen("Genome_C_Imputed", "C", "5860I73K8")


# Filter the genomes to only include common SNPs
counting_again = 0
with open("tgp_plink.map") as infile, open("tgp_samples.map",'w') as outfile:
	onsnps = 0
	offsnps = 0
	for each in infile:
		inline = each.strip().split()
		rsid = inline[1]
		value = int(inline[3])
		
		# If the RSID is not in our map of imputed 23andMe SNPs, remove it
		if rsid not in gen_count and value > 0:
			value = value * -1
			outfile.write(inline[0] + "\t" + inline[1] + "\t" + inline[2] + "\t" + str(value) + "\n")
		# If it is in the map, but not all three individuals have the SNP, remove it
		elif rsid in gen_count and len(gen_count[rsid]) < 3 and value > 0:
			value = value * -1
			outfile.write(inline[0] + "\t" + inline[1] + "\t" + inline[2] + "\t" + str(value) + "\n")
		# Otherwise, write out the original line
		else:
			outfile.write(each)

		# Count the number of valid and invalid SNPs for processing
		if value > 0:
			onsnps = onsnps + 1
		else:
			offsnps = offsnps + 1
	print(onsnps)
	print(offsnps) 
