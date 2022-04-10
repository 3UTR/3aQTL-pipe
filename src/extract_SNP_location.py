'''
SNP location can be extract directly from processed gt file (embeded in SNP id)
run on python2 environment
usage: python extract_SNP_location.py --genotype_bed /path/to/GEUVADIS.all_chrs.gt012.bed --output /path/to/output/snp_location.txt
'''
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--genotype_bed',type=str,help="provide the transformed genotype matrix")
parser.add_argument('--output',type=str, default="snp_location.txt",help="specify SNP location file")

args = parser.parse_args()


fh = open(args.genotype_bed,'r')
fho = open(args.output,'w')
header = fh.readline()
print("SNP\tChr\tPos",file=fho)
for line in fh.readlines():
    line = line.strip()
    snp = line.split("\t")[0]
    w = snp.split("_")
    if len(w)>=2:
        chrom,pos = snp.split("_")[0:2]
        print("%s\t%s\t%s" % (snp,w[0],w[1]), file=fho)
    else:
        print("Error:",snp)
fh.close()
fho.close()

