import property_match_test
import shifting_test
import regression_test
import os
import argparse

 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-ql","--query", default = "Query/Enh_cell",help="Relative folder location of query files")
    parser.add_argument("-al","--annot", default = "Annot", help="Relative folder location of annotation files")
    parser.add_argument("-cl","--confound",default = "table.txt.pruned.gz", help = "Confounding factors file")
    parser.add_argument("-rl","--rsid", default = "GWAS",help="Relative folder location of GWAS files")    
    parser.add_argument("-qi","--query_index",help="Index of query file you want to run test with annotation files(0~114)",type = int)
    parser.add_argument("-t", "--test_type",help = "Type of test you want to run(Property, Regression, Shifting)")
    parser.add_argument("-n", "--nperm", default = 1000, help = "Number of permutation, default = 1000", type = int)    
    args = parser.parse_args()
    
    query_list = os.listdir(os.path.relpath(args.query))
    query_file = query_list[args.query_index]    
    
    annot_list = os.listdir(os.path.relpath(args.annot))
    
    if args.test_type == "Property":
        result = property_match_test.run_one_round_property(args.query+'/'+query_file,args.annot,args.confound,args.nperm)
        result_file = open("Property_test.txt",'a')
        result_file.write(str(args.query_index)+" ")
        for annot_index in range(len(result)):
            result_file.write(str(result[annot_index][1]))
            if annot_index != len(result)-1:
                result_file.write(", ")
        result_file.close()
    elif args.test_type == "Regression":
        result = regression_test.run_one_round_regression(args.query+'/'+query_file,args.annot,args.confound,args.nperm)
        result_file = open("Regression_test.txt",'a')
        result_file.write(str(args.query_index)+" ")
        for annot_index in range(len(result)):
            result_file.write(str(result[annot_index][1]))
            if annot_index != len(result)-1:
                result_file.write(", ")
            else:
                result_file.write("\n")
        result_file.close()
    elif args.test_type == "Shifting":
        result = shifting_test.run_one_round_shifting(args.query+'/'+query_file,args.annot,args.confound,args.nperm)
        result_file = open("Shifting_test.txt",'a')
        result_file.write(str(args.query_index)+" ")
        for annot_index in range(len(result)):
            result_file.write(str(result[annot_index][1]))
            if annot_index != len(result)-1:
                result_file.write(", ")
