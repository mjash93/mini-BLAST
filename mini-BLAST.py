def readfile(filename):
	"""
		uses .read to read the entire file, and the splits it into a list of strings where each 
		string is in essence a paragraph
	"""
	return open(filename, 'r').read().split('\n\n')

def compare_left(query, string, q0, s0):
	"""
		Function is designed to scan left from the start of the k-mer, continuing until
		there is a mismatch in the sequence and query. We return the number of similarities in 
		addition to the original match. 

		INPUTS:
			> query -- string representing the query (a set)
			> string -- string representing the sequence / database being scanned for a match to the query
			> q0 -- starting position of the match in the query
			> s0 -- starting position of the match in the detabase
		OUTPUT:
			> L -- length of additional matches to the left of the k-mer
	"""
	L = 0
	i = 1
	KEY = True
	try:
		while (KEY is True):
			if query[q0 - i] == string[s0 - i]:
				L += 1
				i += 1
			else: KEY = False
		return L
	except IndexError:
		return L 

def compare_right(query, string, q0, s0, num_k):
	"""
		Function is designed to scan right from the end of the k-mer, continuing until
		there is a mismatch in the sequence and query. We return the number of similarities in 
		addition to the original match. 

		INPUTS:
			> query -- string representing the query (a set)
			> string -- string representing the sequence / database being scanned for a match to the query
			> q0 -- starting position of the match in the query
			> s0 -- starting position of the match in the detabase
			> num_k -- length of the k-mer match
		OUTPUT:
			> L -- length of additional matches to the right of the k-mer
	"""
	L = 0
	i = 1
	KEY = True
	try:
		while (KEY is True):
			if query[q0:(q0 + num_k + i)][-1] == string[s0:(s0 + num_k + i)][-1]:
				L += 1
				i += 1
			else: KEY = False
		return L
	except IndexError:
		return L

def all_kmer(string, num_k):
	"""
		Function returns all instances of a kmere in list form. Returns a dict of kmere's with 
		the associated value being a list with the starting positions of each mere. 

		INPUTS:
			> string -- sequence / database to sequence for k-mer's
			> num_k -- number denoting the length of the k-mer to be searched
	"""
	kmer = {}
	for i in range(0 , (len(string)-num_k+1)):
		if string[i:(i+num_k)] not in kmer.keys():
			kmer[string[i:(i+num_k)]] = [i]
		else:
			kmer[string[i:(i+num_k)]].append(i)
	return kmer

def match_kmer(query, string_list, num_k, cutoff, file_out):
	"""
		match_kmer scans a query againsta  database and returns matches along with
		corresponding location and score. 

		Input:
			query = string to check for matches against
			string_list = list of strings (sequences). Our database
			num_k = the length of the kmer to be searched for
			cutoff = the HSP cutoff value
			file_out = name of the file to write results too
		Output:
			Writes results to a file -- for each sequence, results stored in a 
			dict where the key is the kmere and the value is a list of tuples where
			the first element represents the starting location of the match and the second
			represents the HSP value.  
	"""
	f_out = open(file_out, 'w')  # open output file
	query_dict = all_kmer(query, num_k) # returns a dict with list values of all matches for the query string
	
	# We begin the search by selecting each element of 'string_list'. We'll finf all HSP's in each 
	# sequence / database that exceed a threshold and report their location as well as their score.   

	for string in string_list:
		HSP_dict ,first_occ_dict = {}, {}
		max_hsp = num_k
		string = string.upper()

		# We iterate through our database, looking at each entry to see whether  k-mer exists

		for i in range(0,(len(string)-num_k+1)):
			L = num_k # if there is a match, will be at least length k

			# If we locate a match between database and query, we begin the BLAST algorithm and 
			# begin to extend the k-mere to see how long a match we can obtain. 

			if string[i:(i+num_k)] in query_dict.keys():
				# we search not just from the 1st instance of the k-mer in the query, but from each
				# instance
				for sv in range(len(query_dict[string[i:(i+num_k)]])):
					q0 = query_dict[string[i:(i+num_k)]][sv]
					ml = compare_left(query, string, q0, i)
					mr = compare_right(query, string, q0, i, num_k)
					HSP_value = L + mr + ml
					# If your HSP value is greater than the threshold cutoff, categorize and save their scores
					# as well as their start/end positons for output
					if HSP_value >= cutoff:
						# we only want to count a HSP occurance once
						if query[(q0-ml):(q0+num_k+mr)] not in first_occ_dict.keys():
							first_occ_dict[query[(q0-ml):(q0+num_k+mr)]] = ((q0-ml,q0+num_k+mr-1),(i-ml,i+num_k+mr-1))

						if query[(q0-ml):(q0+num_k+mr)] not in HSP_dict.keys():
							HSP_dict[query[(q0-ml):(q0+num_k+mr)]] = 1
						else:
							HSP_dict[query[(q0-ml):(q0+num_k+mr)]] += 1

		if bool(HSP_dict) is not False:
			f_out.write('We found a HSP in sequence: \n\n{}\n\n'.format(string))

			for HSP in HSP_dict:
				f_out.write('We found an occurance of {} at {}.\n'.format(HSP, first_occ_dict[HSP]))



	f_out.close()

def main():
	# ask for all relevant inputs to the problem
	sf = 'database.txt'
	qf = 'query.txt'
	file_out = 'BLAST_OUTPUT.txt'
	k_val = 4 # can hardcode for something else you prefer
	cutoff = 5 # can caustomize this as you wish

	# open and read the contents of the input file, and prepare output file
	query = readfile(qf)[0].upper()
	database = readfile(sf)
	
	match_kmer(query, database, k_val, cutoff, file_out)	

if __name__=='__main__':
	main()