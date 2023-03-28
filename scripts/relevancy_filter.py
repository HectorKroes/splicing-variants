import sys

def find_start(lines):
	for i in range(len(lines)):
		if lines[i][0] == "#":
			continue
		else:
			return i

def verify_annotation_presence(line):
	squirls = line.find("SQUIRLS_SCORE=")
	spliceai = line.find("SpliceAI=")
	return squirls, spliceai

def score_formatter(score):
	if score == '.':
		return 0.0
	else:
		return float(score)

def verify_relevance(line, index, spliceai_cutoff, squirls_cutoff):
	annotations = line[index:].split('\t')[0].split('=')
	squirls_score = score_formatter(annotations[1].split(';')[0])
	max_spliceai_score = max([score_formatter(i) for i in annotations[2].split('|')[2:6]])
	if squirls_score >= squirls_cutoff and max_spliceai_score >= spliceai_cutoff:
		return True
	else: 
		return False

def save_results(results, file_name):
	with open(file_name, 'w') as file:
		file.write("".join(results))

def filter_variants(variants, start, results, spliceai_cutoff, squirls_cutoff, file_name):
	for i in range(len(variants)):
		score1, score2 = verify_annotation_presence(variants[i])
		if score1 != -1 and score2 != -1:
			if verify_relevance(variants[i], score1, spliceai_cutoff, squirls_cutoff):
				results.append(variants[i])
	save_results(results, file_name)

if __name__ == "__main__":
    file_path = sys.argv[1]
    spliceai_cutoff = float(sys.argv[2])
    squirls_cutoff = float(sys.argv[3])
    results_file = f"relevant_{''.join(file_path.split('_', 1)[1])}"
    with open(file_path) as file:
    	lines = file.readlines()
    start = find_start(lines)
    header = lines[:start]
    variants = lines[start:]
    filter_variants(variants, start, header, spliceai_cutoff, squirls_cutoff, results_file)