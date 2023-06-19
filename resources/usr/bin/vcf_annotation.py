import sys

def find_start(lines):
	for i in range(len(lines)):
		if lines[i][0] == "#":
			continue
		else:
			return i

def header_formatter(header):
	for i in range(len(header)):
		if header[i] == '##FILTER=<ID=SQUIRLS,Description="Squirls considers the variant as pathogenic if the filter is present">\n':
			header[i] = '##INFO=<ID=SplicePrediction,Number=A,Type=String,Description="Splicing pathogenicity prediction">\n'
	return header

def verify_annotation_presence(line):
	squirls = line.find("SQUIRLS_SCORE=")
	spliceai = line.find("SpliceAI=")
	return squirls, spliceai

def score_formatter(score):
	if score == '.':
		return 0.0
	else:
		return float(score)

def verify_relevance_both(line, index, spliceai_cutoff, squirls_cutoff, mode):
	annotations = line[index:].split('\t')[0].split('=')
	squirls_score = score_formatter(annotations[1].split(';')[0])
	max_spliceai_score = max([score_formatter(i) for i in annotations[2].split('|')[2:6]])
	if mode == "AND":
		if squirls_score > squirls_cutoff and max_spliceai_score > spliceai_cutoff:
			return 'P'
		else: 
			return 'N'
	elif mode == 'OR':
		if squirls_score > squirls_cutoff:
			return 'P'
		elif max_spliceai_score > spliceai_cutoff:
			return 'P'
		else:
			return 'N'

def verify_relevance_squirls(line, index, squirls_cutoff):
	annotations = line[index:].split('\t')[0].split('=')
	squirls_score = score_formatter(annotations[1].split(';')[0])
	if squirls_score > squirls_cutoff:
		return 'P'
	else:
		return 'N'

def save_results(results, file_name):
	with open(file_name, 'w') as file:
		file.write("".join(results))

def annotate(variant_line, prediction):
	columns = variant_line.split('\t')
	columns[6] = '.'
	columns[7] = columns[7].replace('\n', '')
	if len(columns) == 8:
		columns[7] = f'{columns[7]};SplicePrediction={prediction}\n'
	else:
		columns[7] = f'{columns[7]};SplicePrediction={prediction}'
	return '\t'.join(columns)

def filter_variants(variants, start, results, spliceai_cutoff, squirls_cutoff, file_name, mode):
	for i in range(len(variants)):
		score1, score2 = verify_annotation_presence(variants[i])
		if score1 != -1 and score2 != -1:
			results.append(annotate(variants[i], verify_relevance_both(variants[i], score1, spliceai_cutoff, squirls_cutoff, mode)))
		elif score1 != -1:
			if mode == "OR":
				results.append(annotate(variants[i], verify_relevance_squirls(variants[i], score1, squirls_cutoff)))
			elif mode == "AND":
				results.append(variants[i])
		else:
			results.append(variants[i])
	save_results(results, file_name)

if __name__ == "__main__":
    file_path = sys.argv[1]
    spliceai_cutoff = float(sys.argv[2])
    squirls_cutoff = float(sys.argv[3])
    mode = str(sys.argv[4]).upper()
    results_file = f"splice_{file_path.split('/')[-1]}"
    with open(file_path) as file:
    	lines = file.readlines()
    start = find_start(lines)
    header = header_formatter(lines[:start])
    variants = lines[start:]
    filter_variants(variants, start, header, spliceai_cutoff, squirls_cutoff, results_file, mode)