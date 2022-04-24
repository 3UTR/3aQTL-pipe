import argparse
import os.path

# -- Functions
def load_sample_list(input_sample_list_file):
    sample_list = []
    for line in open(input_sample_list_file,'r'):
        line = line.strip()
        sample_id = line.split("\t")[0]
        sample_list.append(sample_id)

    return sample_list

def extract_total_reads(input_flagstat_file):
    num_line = 0
    total_reads = '-1'
    #print input_flagstat_file
    for line in open(input_flagstat_file,'r'):
        num_line += 1
        if num_line == 5:
            total_reads = line.strip().split(' ')[0]
            break
    return total_reads


# -- Main

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sample_list',help="the file contains all samples")
    parser.add_argument('--path_flagstat',default="./tmp",help="the location of flagstat files")
    parser.add_argument('--path_wig',default="./wig",help="the location of wig files")
    parser.add_argument('--output',help="the final output file with read depth of each sample")

    args = parser.parse_args()

    selected_samples = load_sample_list(os.path.abspath(args.sample_list))
    path_wig = os.path.abspath(args.path_wig)
    path_flagstat = os.path.abspath(args.path_flagstat)
    
    if path_wig[-1] != "/":
        path_wig += "/"
    else:
        pass

    if path_flagstat[-1] != "/":
        path_flagstat += "/"
    else:
        pass

    fho = open(args.output,'w')
    for sample in selected_samples:
        filename_sample = sample + ".flagstat"
        read_depth = extract_total_reads(path_flagstat + filename_sample)
        wig_location = path_wig + sample + ".wig"
        print("%s\t%s" % (wig_location,read_depth),file=fho)

    fho.close()
