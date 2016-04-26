# Personal Genome Script 1A: samples_to_hapmap.py 

# Helper function to synchronize SNPs
def filter_missing_snps(snpSrc, snpRef):
	to_remove = []
	for each in snpSrc:
		if each not in snpRef:
			to_remove.append(each)
	for each in to_remove:
		del snpSrc[each]
	return snpSrc


# Helper Function to process a sample genome file
def process_sample(filename):
	snpMap = {}
	with open(filename) as sample_file:
		for each in sample_file:
			# If not a comment
			if each[0] != "#":
				# Parse the line and map rsid -> genotype call
				inline = each.strip().split()
				snpMap[inline[0]] = inline[3]
	return snpMap


# Helper Function to write a sample to a .ped file
def write_to_ped(pedfile, sample_id, sample_reads):
	pedfile.write(sample_id + "\t" + sample_id + "\t0\t0\t1\t-9\t")
	# For each SNP in the .map file
	with open("samples_hapmap.map") as mapfile:	
		for each in mapfile:
			rsid = each.strip().split()[1]
			# If the rsid is in the sample and has a call
			if rsid in sample_reads and '-' not in sample_reads[rsid]:
				call = sample_reads[rsid]
				# If the call has both alleles, write them to the .ped file
				if len(call) > 1:
					pedfile.write(call[0] + " " + call[1] + "\t")
				# If the call has only one allele, write as if homozygous
				else:
					pedfile.write(call[0] + " " + call[0] + "\t")
			# If the rsid is not in the sample, then write 0 0 to indicate a missing call
			else:
				pedfile.write("0 0\t")
		pedfile.write("\n")



# Process sample genomes
snpA = process_sample("genome_A.txt")
print(len(snpA))
snpB = process_sample("genome_B.txt")
print(len(snpB))
snpC = process_sample("genome_C.txt")
print(len(snpC))


# Filter sample maps to include only shared SNPs
snpA = filter_missing_snps(snpA, snpB)
snpA = filter_missing_snps(snpA, snpC)
snpB = filter_missing_snps(snpB, snpA)
snpC = filter_missing_snps(snpC, snpA)


rscheck = {}
with open("hapmap3.map") as mapfile:
	count = 0
	for each in mapfile:
		mapline = each.strip().split()
		# If the RSID is in our shared SNP set
		if mapline[1] in snpA:
			# Map RSID -> .map position
			rscheck[mapline[1]] = count
		count = count + 1
print("rs check length:")
print(len(rscheck))


# Remove SNPs not present in the hapmap file
filter_missing_snps(snpA, rscheck)
filter_missing_snps(snpB, rscheck)
filter_missing_snps(snpC, rscheck)


# Write new map file, filtering out non-common SNPs
print("writing samples_hapmap")
hapmap_lines = []
with open("hapmap3.map") as infile, open("samples_hapmap.map", 'w') as outfile:
	onsnps = 0
	offsnps = 0
	for each in infile:
		inline = each.strip().split()
		value = int(inline[3])
		# If the SNP is a common SNP, then write back the original line
		if inline[1] in rscheck:
			outfile.write(each)
		# If the SNP is not a common SNP, replace the base-position field with a negative number
		else:
			value = value * -1	
			newline = inline[0] + "\t" + inline[1] + "\t" + inline[2] + "\t" + str(value) + "\n"
			outfile.write(newline)

		# Count the number of valid and invalid SNPs for processing
		if value > 0:
			onsnps = onsnps + 1
		else:
			offsnps = offsnps + 1
	print(onsnps)
	print(offsnps)				


# Write the sample genome ped file
print("writing ped")
with open("samples_hapmap.ped", 'w') as pedfile:
	write_to_ped(pedfile, "A", snpA)
	write_to_ped(pedfile, "B", snpB)
	write_to_ped(pedfile, "C", snpC)
