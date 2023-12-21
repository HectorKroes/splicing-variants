import sys

def find_start(lines):
	for i in range(len(lines)):
		if lines[i][0] == "#":
			continue
		else:
			return i

def save_results(results, file_name):
	with open(file_name, 'w') as file:
		file.write("".join(results))

def filter_bool(variant_line):
	columns = variant_line.split('\t')
	info = columns[7]
	gnomad_af = float(info.split('gnomad_popmax_af=')[1].split(';')[0])
	if gnomad_af < float(af_filter):
		return True
	else:
		return False

if __name__ == "__main__":

    file_path = sys.argv[1]
    do_filter = sys.argv[2]
    af_filter = sys.argv[3]

    if do_filter == 'true':
	    with open(file_path) as file:
	    	lines = file.readlines()
	    start = find_start(lines)
	    header = lines[:start]
	    variants = lines[start:]
	    passing_indices = []
	    index = 0
	    for variant in variants: 
	    	if filter_bool(variant):
	    		passing_indices.append(index)
	    	index += 1
	    for passing_index in passing_indices:
	    	header.append(variants[passing_index])
	    save_results(header, 'slivar_vcf.vcf')
