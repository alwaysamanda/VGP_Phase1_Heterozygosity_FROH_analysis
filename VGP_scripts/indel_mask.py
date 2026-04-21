import numpy 
import matplotlib.pyplot as plt
import gzip
import re
import argparse

parser=argparse.ArgumentParser(description="Input Genotype Likelihood VCF and output BED")

parser.add_argument('--input_gen_vcf',type=str,required=True,help="Input file path")
parser.add_argument('--Chromosome_name',type=str,required=True,help="Name of the chromosome")
parser.add_argument('--output_BED',type=str,required=True,help="Ouptut Path to store multihet")

args=parser.parse_args()

gvcf_file=args.input_gen_vcf
chromosome_name=args.Chromosome_name
output_bed=args.output_BED

positions=[]
with gzip.open(gvcf_file,"rt") as file:
    for line in file:
        if (line[0]=="#"):
            continue  
        current_line=line.strip().split('\t')
        pos=int(current_line[1])
        info=current_line[7]
        dp_match = re.search(r"DP=(\d+)", info)
        if dp_match and int(dp_match.group(1)) == 0:
            positions.append(pos)

start_window=0
k=0

with open(output_bed, "w") as output_bed_file:
    for i in range(0, len(positions) - 1):
        if int(positions[i]) + 1 == int(positions[i + 1]):
            if k == 0:
                start_window = positions[i]
            k += 1
        else:
            if k == 0:
                output_bed_file.write(f"{chromosome_name}\t{positions[i]}\t{positions[i]}\n")
            else:
                output_bed_file.write(f"{chromosome_name}\t{start_window}\t{positions[i]}\n")
                k = 0  





