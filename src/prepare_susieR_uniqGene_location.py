'''
Desc:
Obtain unique and significant aGenes and their extended location (specified by --extend_size, default:1Mb)
Input files: 3utr_location.txt, Cis_aQTL_all_control_gene_exprs.txt
Output format: Gene\tchr:start-end
'''

import argparse
import os.path
import time

# - Main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--utr_loc_file',help="input the reference 3'UTR location file, e.g. 3utr_location.txt")
    parser.add_argument('--aQTL_map',help="specify the aQTL mapping file")
    parser.add_argument('--Max_FDR',type=float,default=0.05,help="specify the maximum of FDR to select significant aQTL associations")
    parser.add_argument('--extend_size',type=int, default=1000000,help="Int, extend N bp (default N=1e6) at both sides")
    parser.add_argument('--outdir',help="specify the output dir")
    parser.add_argument('--output',help="specify the output file name")

    args = parser.parse_args()

    # -- extract gene and location
    print("Start extracting 3'UTR location file...")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    fh = open(args.utr_loc_file,'r')
    gene2loc = {}
    head = fh.readline()
    ext_size = int(args.extend_size)

    for line in fh.readlines():
        line = line.strip()
        w = line.split("\t")
        start = max([0,int(w[2]) - ext_size])
        end = int(w[3]) + ext_size
        loc = w[1] + ":" + str(start) + "-" + str(end)
        gene2loc[w[0]] = loc
    fh.close()

    # -- processing aQTL mapping files
    print("Start processing aQTL mapping file...")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    fdr_cutoff = float(args.Max_FDR)
    fh = open(args.aQTL_map,'r')
    outdict = {}
    for line in fh.readlines()[1:]:
        line = line.strip()
        w = line.split("\t")
        gene = w[1]
        fdr = float(w[5])
        if fdr < fdr_cutoff:
            if gene not in outdict:
                outdict[gene] = gene2loc[gene]
            else:
                continue
        else:
            continue

    fh.close()

    print(len(outdict),"unique aGenes are processed.")
    

    output_path = os.path.abspath(args.outdir)
    fho = open(output_path + "/" + args.output,'w')
    i = 0
    for gene in outdict:
        i += 1
        print(gene + "\t" + outdict[gene], file=fho)
        print(i)
    fho.close()
    print("Done!")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
