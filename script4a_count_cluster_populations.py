# Personal Genome Script 4A: count_cluster_populations.py

import sys, os

id_map = {}
with open("relationships_w_pops_121708.txt") as infile:
	lines = infile.readlines()
	for each in lines:
		line = each.strip().split()
		full_id = line[0] + "_" + line[1]
		id_map[full_id] = line[6]
	infile.close()

cluster_map = {}
with open("plink.cluster1") as infile:
	lines = infile.readlines()
	for each in lines:
		line = each.strip().split()
		cluster_map[line[0]] = []
		for person in line[1:]:
			if person not in id_map:
				print(person)
			cluster_map[line[0]] = cluster_map[line[0]] + [person] 
	infile.close()


pop_map = {}
count = 0
with open("cluster_count.txt", 'w') as outfile:
	for each in cluster_map:
		id_list = cluster_map[each]
		pop_count = {}
		
		for person in id_list:
			pop_name = "Unknown"
			if person in id_map:
				pop_name = id_map[person]
			else:
				count = count + 1
			if pop_name not in pop_count:
				pop_count[pop_name] = 0
			pop_count[pop_name] = pop_count[pop_name] + 1
			pop_map[each] = pop_count
	
	sorted_keys = list(pop_map.keys())
	sorted_keys.sort(key = lambda a: int(a.split("-")[1]))	
	for each in sorted_keys: 
		outfile.write(each + ":\n")
		for e in pop_map[each]:
			outfile.write(e + ":\t" + str(pop_map[each][e]) + ", ")
		outfile.write("\n")
	outfile.close()

print(count)
