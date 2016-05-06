# Personal Genome Script 3A: retrieve_target_samples.py
# Requires target_pop.txt file

# Retrieve the list of target CEU and TSI samples
pop_dict = {}
with open("target_pop.txt") as popfile:
	for each in popfile:
		popline = each.strip().split()
		name = popline[0] + "_" + popline[1]
		pop_dict[name] = popline[6]

# Hapmap .map file dictionary: maps the rsid index -> rsid
hapmap_dict = {}
print("hapmap processing")
with open("hapmapandme.map") as mapfile:
	count = 0
	for each in mapfile:
		mapline = each.strip().split()
		hapmap_dict[count] = mapline[1]	
		count = count + 1
num_genes = len(hapmap_dict)

# TGP .map file dictionary: maps the rsid index -> rsid
tgp_dict = {}
print("tgp map processing")
with open("tgp_plink.map") as infile:
	count = 0
	for each in infile:
		tgp_line = each.strip().split()
		tgp_dict[count] = tgp_line[1]
		count = count + 1
num_tgp_genes = len(tgp_dict)

# Open the original .ped file and write to a new .ped file
print("ped processing")
ind_count = 0
with open("hapmapandme.ped") as pedfile, open("combined_pop.ped", 'w') as outfile:
	for each in pedfile:
		pedline = each.strip().split()
		name = pedline[0] + "_" + pedline[1]
	
		# If this individual is in our target sample set
		if name in pop_dict:
			ind_count = ind_count + 1
			
			# Create a new dictionary for the individual's genomes: rsid -> genotype call
			new_person_genome = {}
			for i in range(num_genes):
				new_person_genome[hapmap_dict[i]] = pedline[(i*2)+6] + pedline[(i*2)+7]

			# Write a new line to the output file in TGP .map format	
			famname = pedline[0]
			idname = pedline[1]
			outfile.write(famname + "\t" + idname + "\t0\t0\t0\t-9\t")

			valid_count = 0
			# For each SNP in the tgp genome dictionary
			for i in range(num_tgp_genes):
				rsid = tgp_dict[i]
				# If the SNP is in the individual's genome and the call has two alleles
				if rsid in new_person_genome and len(new_person_genome[rsid]) == 2:
					# Write the SNP call to the new .ped file
					snp_call = new_person_genome[rsid]
					outfile.write(snp_call[0] + " " + snp_call[1] + "\t")
					valid_count = valid_count + 1
				# Else, write 0 0 for unknown genotype
				else:
					outfile.write("0 0\t")
			print(valid_count)
			outfile.write("\n")
