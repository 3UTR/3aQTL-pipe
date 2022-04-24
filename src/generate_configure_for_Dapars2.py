import argparse
import os.path

# -- Functions
def extract_all_wigs(input_seq_depth_file):
    wig_files = []
    for line in open(input_seq_depth_file,'r'):
        line = line.strip()
        w = line.split("\t")
        wig_files.append(w[0])

    return wig_files


# -- Main

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--annotation_3utr',help="the location of reference 3'UTR bed file")
    parser.add_argument('--wigFile_depth',help="the index file contains all wig files and read depth")
    parser.add_argument('--coverage_threshold',type=str,default="10",help="specify the threshold of coverage,default=10")
    parser.add_argument('--threads',type=str,default="1",help="specify the number of threads used,default=1")
    parser.add_argument('--out_dir_prefix',type=str, default="Dapars2_out",help="specify the directory prefix of dapars2 output,default='Dapars2_out'")
    parser.add_argument('--out_file_prefix',type=str,default="Dapars2",help="specify the name of result file prefix of dapars2,default='Dapars2'")
    parser.add_argument('--out_config_name',type=str,default="Dapars2_running_configure.txt",help="specify configure file name,default='Dapars2_running_configure.txt'")

    args = parser.parse_args()

    configure_file_name = args.out_config_name
    fho = open(configure_file_name,'w')
    # print Annotated_3UTR
    print("# Specify the reference of 3'UTR region", file=fho)
    print("\nAnnotated_3UTR=" + os.path.abspath(args.annotation_3utr), file=fho)

    # print wig files
    all_wig_files = extract_all_wigs(os.path.abspath(args.wigFile_depth))
    print("\n# A comma separated list of wig files of all samples", file=fho)
    print("\nAligned_Wig_files=" + ",".join(all_wig_files), file=fho)
    
    # specify Output_directory and Output_result_file
    print("\nOutput_directory=" + args.out_dir_prefix, file=fho)
    print("\nOutput_result_file=" + args.out_file_prefix, file=fho)

    # specify Coverage_threshold
    print("\n# Specify Coverage threshold", file=fho)
    print("\nCoverage_threshold=" + args.coverage_threshold, file=fho)

    # specify the Num_Threads to process the analysis
    print("\n# Specify the number of threads to process the analysis", file=fho)
    print("\nNum_Threads=" + args.threads, file=fho)

    # Provide sequencing_depth_file for normalization
    print("\n# Provide sequencing depth file for normalization", file=fho)
    print("\nsequencing_depth_file=" + os.path.abspath(args.wigFile_depth), file=fho)
    fho.close()
    

