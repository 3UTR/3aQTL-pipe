'''
after getting *.frq and *.FORMAT files by vcftools,
this script will extracting allelic gt info from *.frq and recode gt code in *.FORMAT into 012 format
'''

import argparse
import time


# - Functions
# extracting allelic genotype from frq
def extract_gt_from_frq(frq_file):
    snp2gt = {}
    fh = open(frq_file,'r')
    for line in fh.readlines()[1:]:
        line = line.strip()
        w = line.split("\t")
        if "chr" not in w[0]:
            snp = "chr" + w[0] + "_" + w[1]
        else:
            snp = w[0] + "_" + w[1]
        allele_ref = w[4].split(":")[0]
        allele_alt = w[5].split(":")[0]

        if snp not in snp2gt:
            snp2gt[snp] = allele_ref + "_" + allele_alt
        else:
            continue
    fh.close()
    return snp2gt



# recode genotype to 012 code
def recode_with_012(gt):
    if "." in gt:
        return "NA"
    else:
        allele_1 = int(gt[0])
        allele_2 = int(gt[-1])
        return str(allele_1+allele_2)

# - Main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--frq',help="input the frq file generated by vcftools")
    parser.add_argument('--GT',help="specify the GT format file generated by vcftools")
    parser.add_argument('--output',help="specify the output file with 012 recoded genotype")

    args = parser.parse_args()

    # -- extract allelic gt from frq file
    print("Start processing frq file...")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    snv_gt = extract_gt_from_frq(args.frq)
    print("Obtain genotype of %d SNPs" % (len(snv_gt)))

    # -- recode gt in GT file
    print("Start processing GT file...")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    fh = open(args.GT,'r')
    fho = open(args.output,'w')
    header = fh.readline().strip().split("\t")
    print("%s\t%s" % ("id","\t".join(header[2:])), file=fho)
    for line in fh.readlines():
        line = line.strip()
        w = line.split("\t")
        if "chr" not in w[0]:
            snp = "chr" + w[0] + "_" + w[1]
        else:
            snp = w[0] + "_" + w[1]
        
        if snp in snv_gt:
            snp = snp + "_" + snv_gt[snp]
        else:
            snp = snp

        gt_012 = list(map(recode_with_012,w[2:]))
        
        print("%s\t%s" % (snp,"\t".join(gt_012)), file=fho)

    fh.close()
    fho.close()
    
    print("recoded gt file has been write to",args.output)
    print("Done!")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
