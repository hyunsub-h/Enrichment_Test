import gzip
import os
import math

def merge_interval_dict(query):
# read query and make an dictionary of unoverlapped interval    
# query, annot: chromosome, starting, ending
# dict: key: chromosome value: sorted list of starting, end(inclusive) as list (starting position, True) (ending position, False)
    interval_dict = {}
    for item in query:
        interval_dict.setdefault(item[0],[])
        interval_dict[item[0]]= interval_dict[item[0]] + [(item[1],True),(item[2],False)]
    for chrom in interval_dict:
        interval_dict[chrom]= sorted(interval_dict[chrom],key = lambda x:x[0])
    for chrom in interval_dict:
        position_list = interval_dict[chrom]
        starting_positions = [] # unoverlapped interval starting position in order
        ending_positions = [] # unoverlapped interval ending position in order
        counting_value = 0 # # for ending - # of starting so far
        for item in position_list:
            if item[1]: # if starting position
                if counting_value == 0:
                    starting_positions.append(item[0])
                counting_value -=1
            else:
                counting_value += 1
                if counting_value == 0:
                    ending_positions.append(item[0])
        interval_dict[chrom] = [starting_positions,ending_positions]
    return interval_dict
            
def make_confound_dict(confound_txt):
# key: rsid, value: list of [chromosome, position, confounding value to compare]
    confound_dict = {}
    with gzip.open(confound_txt) as confound_file:
        for line in confound_file:
            tmp = line.strip().split()
            tmp = tmp[0].split('|')[:-1]+tmp[1:]
            confound_value_to_compare = 0
            multiplier = 1
            for confound_value in  tmp[4:]:
                confound_value_to_compare += float(confound_value)*10**multiplier
                multiplier+=1
            confound_dict[tmp[0]] = [tmp[1],int(tmp[2]),confound_value_to_compare] 
    return confound_dict

def make_rsid_info(rsid_folder):
# key:rsid, value: dictionary of {key:disease value:-log(p_value)}
    file_list = os.listdir(os.path.relpath(rsid_folder))
    rsid_dict = {}
    with open(os.path.relpath(rsid_folder+'/background.rsid.lp')) as background_rsid:
        for line in background_rsid:
            rsid_dict.setdefault(line.strip(),{})
    for filename in file_list:
        if filename != 'background.rsid.lp':
            with open(os.path.relpath(rsid_folder+'/'+filename)) as disease_p:
                for line in disease_p:
                    p_info = line.strip().split()
                    rsid_dict.setdefault(p_info[0],{})
                    rsid_dict[p_info[0]][filename[:-8]]=-math.log(float(p_info[1])) if p_info[1]!='0' else 0#float('inf')
    return rsid_dict 
