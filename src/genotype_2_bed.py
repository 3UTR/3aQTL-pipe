# Transform the 012 genotype matrix into a bed file
# by adding three columns (chr start end) as the first three columns
# output the Header line of genotype matrix for further use

import argparse
# - Main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--genotype', help="provide the genotype matrix used in 3'aQTL mapping")
    parser.add_argument('--out_bed', help="specify the output bed file name")
    parser.add_argument('--out_header', help="specify a file to store the genotype header")

    args = parser.parse_args()

    fh = open(args.genotype,'r')
    fho = open(args.out_bed,'w')
    fho_header = open(args.out_header,'w')
    header = fh.readline().strip()
    print(header, file=fho_header)
    fho_header.close()
    i = 0
    for line in fh.readlines():
        i += 1
        line = line.strip()
        snp = line.split("\t")[0]
        chrom,pos,ref,alt,release = snp.split("_")
        pos = int(pos)
        print("%s\t%d\t%d\t%s" % (chrom,pos-1,pos,line), file=fho)
#        print(i)

    print(i,"SNPs been processed!")
fh.close()
fho.close()
