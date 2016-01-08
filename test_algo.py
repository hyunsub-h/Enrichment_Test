import random

# format of the file 
# BED FILE : Name, starting position, ending position, name of BED line, score, strand, thickstart, thickend, itemrgb, blockcount, blocksizes, blockstarts
# GRange : seqname, ranges, strand, metadata
# ranges : represented as (a,b) where its range is range(a,b) = (a,a+1,...,b-1)
# For this program, we will work with GRange format with score as the first metadata
# i.e. for item in query(or annot), item[3]:score


#
# property matching method
#

property_test(query, annot, confounding_factors){
	
}

#
# regression method
#

regression_test(qeury, annot, confounding_factors)

#
# shifting method
#

# version 3 With locus & circulization


make_annot_dict(annot){
	annot_dict = {}
	for item in annot_dict:
		annot_dict[item[0]]=query_dict.get(item[0],[])+[('s',item[1][0]),('e',item[1][-1]+1)] #key: seq_num value:starting position and ending position(inclusive) of ranges of GRange objects whose seq_name is key
	for item in annot_dict:
		annot_dict[item].sort(key = lambda x : x[1])
	return annot_dict
}


count_query_overlap_with_dict(query, annot_dict){
	query_dict = {}
	for item in query:
		query_dict[item[0]]=query_dict.get(item[0],[])+[item[1][0]] #key: seq_num value:starting position of ranges of GRange objects whose seq_name is key
	count = 0 
	for item in query_dict:
		query_pos_s = query_dict[item]
		query_pos_s.sort()
		annot_pos = annot_dict[item]
		annot_count = 0 #count # of end - # of start
		query_counter_pos = 0
		annot_counter_pos = 0
		while query_counter_pos < len(query_pos_s):
			while annot_pos[annnot_counter_pos][1]<query_pos_s[query_counter_pos]:
				if annot_pos[annot_counter_pos][0]=='s':
					annot_count += 1
				else:
					annot_count -= 1
				annot_counter_pos+=1
			if annot_count>0:
				count+=1
			query_counter_pos += 1
	return count

}

cal_annot_range(annot_dict){
	starting_pos = []
	ending_pos = []
	for item in annot_dict:
		starting_pos.append(annot_dict[item][0][1])
		ending_pos.append(annot_dict[item][-1][1])
	starting_pos.sort()
	ending_pos.sort()
	return (starting_pos[0],ending_pos[-1])
}

shifting_test(query, annot, nperm=1000, expand=2){
	random.seed()
	annot_dict = make_annot_dict(annot)
	annot_range = cal_annot_range(annot_dict)
	annot_size = annot_range[1]-annot_range[0]+1 #ending pos was inclusivev
	test_range = (max(0,annot_range[0]-expand * annot_size/2),annot_range[1]+expand * annot_size/2) #inclusive
	observed_count = count_query_overlap_with_dict(query,annot_dict)
	perm_count = []
	for i in range(nperm):
		delta_size = expand*2*(random.random()-0.5)*annot_size
		rand_query = []
		for item in rand_query:
			if item[1][0]+delta_size < test_range[0] and item[1][0] > test_range[0]:
				rand_query.append([item[0],(test_range[1]-(test_range[0]-item[1][0]-delta_size)+1, test_range[1])])
				rand_query.append([item[0],(test_range[0], item[1][1]+delta_size)])
			elif item[1][1]+delta_size > test_range[1] and item[1][1] < test_range[1]:
				rand_query.append([item[0],(item[1][0]+delta_size,test_range[1])])
				rand_query.append([item[0],(test_range[0], test_range[0]+(-test_range[1]+item[1][1]+delta_size)-1)])
			else:
				rand_query.append([item[0],(item[1][0]+delta_size, item[1][1]+delta_size)])
		perm_count.append(count_query_overlap_with_dict(rand_query,annot_dict))
	count = 0
	for num in perm_count:
		if num > observed_count:
			count+=1
	return (observed.count,count/float(nperm))
}

