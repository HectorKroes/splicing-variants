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

def normalize(variant_line):
	columns = variant_line.split('\t')
	columns[3] = columns[3].upper()
	columns[4] = columns[4].upper()
	return '\t'.join(columns)

if __name__ == "__main__":
    file_path = sys.argv[1]
    with open(file_path) as file:
    	lines = file.readlines()
    start = find_start(lines)
    header = lines[:start]
    variants = lines[start:]
    for variant in variants: 
    	header.append(normalize(variant))
    save_results(header, 'normalized_vcf.vcf')
