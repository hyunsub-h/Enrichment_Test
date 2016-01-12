import gzip
import os
import numpy
import bisect
import utils

    
#
# property matching method
#        
    
def property_test(query, annot, confound_dict, nperm): 
# query, annot: chromosome, starting, ending
# confound_dict: key: rsid, value: list of [chromosome, position, confounding value to compare]
# for snp in query -> make another snp with the same confounding value and see if there is matching
    numpy.random.seed()
    query_dict = utils.merge_interval_dict(query)
    annot_dict = utils.merge_interval_dict(annot)
    total_snp = {} # chr: {confounding_value:count}
    query_snp = {} # chr: {confounding_value:count}
    annot_snp = {} # chr: {confounding_value:count} 
    observed_count = 0
    for rsid in confound_dict:
        in_query = 0 # is rsid in query?
        in_annot = 0 # is rsid in annot?
        info = confound_dict[rsid]
        total_chrom_dict = total_snp.setdefault(info[0],{})
        total_chrom_dict[info[2]]=total_chrom_dict.get(info[2],0)+1
        query_pos = query_dict[info[0]]
        annot_pos = annot_dict[info[0]]
        query_index = bisect.bisect(query_pos[0],info[1])-1
        if query_index != -1 and query_pos[1][query_index]>=info[1]:
            query_chrom_dict = query_snp.setdefault(info[0],{})
            query_chrom_dict[info[2]]=query_chrom_dict.get(info[2],0)+1
            in_query = 1
        annot_index = bisect.bisect(annot_pos[0],info[1])-1
        if annot_index != -1 and annot_pos[1][annot_index]>=info[1]:            
            annot_chrom_dict = annot_snp.setdefault(info[0],{})
            annot_chrom_dict[info[2]]=annot_chrom_dict.get(info[2],0)+1
            in_annot = 1
        if in_query*in_annot == 1:
            observed_count +=1
    permutation_count = [0]*nperm
    for chrom in query_snp:
        for confound_value in query_snp[chrom]:
            count_for_given_confound = numpy.random.hypergeometric(annot_snp.get(chrom,{}).get(confound_value,0),total_snp[chrom][confound_value]-annot_snp.get(chrom,{}).get(confound_value,0),query_snp[chrom][confound_value],nperm)
            permutation_count = [sum(x) for x in zip(permutation_count, count_for_given_confound)]
    p_value = 0
    for count in permutation_count:
        if count>observed_count:
            p_value+=1
    p_value /= float(nperm)
    return (observed_count,p_value)
        
#def run_property_test(query_file_name,annot_file_name,confound_txt, nperm=1000):
#    query = []
#    with gzip.open(query_file_name,'rb') as f_query:
#        for line in f_query:
#            to_add = line.strip().split()
#            query.append([to_add[0],int(to_add[1]),int(to_add[2])])
#    annot = []
#    with gzip.open(annot_file_name,'rb') as f_annot:
#        for line in f_annot:            
#            to_add = line.strip().split()
#            annot.append([to_add[0],int(to_add[1]),int(to_add[2])])
#    confound_dict = make_confound_dict(confound_txt)
#    return property_test(query,annot,confound_dict,nperm)

def run_one_round_property(query_file,annot_path,confound_txt,nperm=1000):
    result = []
    confound_dict = utils.make_confound_dict(confound_txt)
    query = []
    with gzip.open(query_file,'rb') as f_query:
        for line in f_query:
            to_add = line.strip().split()
            query.append([to_add[0],int(to_add[1]),int(to_add[2])])
    file_list = os.listdir(os.path.relpath(annot_path))
    for filename in file_list:
        with gzip.open(os.path.relpath(annot_path+'/'+filename)) as f_annot:
            annot = []
            for line in f_annot:            
                to_add = line.strip().split()
                annot.append([to_add[0],int(to_add[1]),int(to_add[2])])
            result.append(property_test(query,annot,confound_dict,nperm))
    return result
