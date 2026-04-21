import msprime
import matplotlib.pyplot as plt 
# import seaborn as sns
import statistics
import numpy as np
import argparse

parser=argparse.ArgumentParser(description="Input file, chromosome name and Output file")

parser.add_argument('--input_file',type=str,required=True,help="Input file path")
parser.add_argument('--Chromosome_name',type=str,required=True,help="Name of the chromosome")
parser.add_argument('--output_BED',type=str,required=True,help="Ouptut Path to store multihet")

args=parser.parse_args()

file=args.input_file
chrom=args.Chromosome_name
output_bed_directory=args.output_BED

counter=0
het=[]
with open(file,"r") as file_reader:
    file_content=file_reader.readlines()
for i,line in enumerate(file_content):
        line=line.split("\t")
        if i==0:
           start_position=int(line[1])
        position=int(line[1])
        if position > (start_position+10000):
            het.append(counter)
            difference=position-start_position
            num_of_window=difference//10000
            if num_of_window!=0:
                for _ in range(num_of_window-1):
                    het.append(0)
            counter = 1
            start_position=start_position+(10000*num_of_window)
            continue
        counter=counter+1

het.append(counter)

percentile_cutoff=95
percentile_list=np.percentile(het,95)

print(f"The 95 percentile Depth for the cutoff : {percentile_list} \n" , flush=True)


start_position_mask=[]
end_position_mask=[]
counter=0
for i,line in enumerate(file_content):
        line=line.split("\t")
        if i==0:
           start_position=int(line[1])
           start_position_mask.append(0)
           end_position_mask.append(start_position)
        position=int(line[1])
        if position > (start_position+10000):
            if counter > percentile_list:
                start_position_mask.append(start_position)
                end_position_mask.append(start_position+10000)
            difference=position-start_position
            num_of_window=difference//10000
            counter = 1
            start_position=start_position+(10000*num_of_window)
            continue
        counter=counter+1

with open(f"{output_bed_directory}/{chrom}_{percentile_cutoff}.txt","w") as bed_writer:
        for start, end in zip(start_position_mask, end_position_mask):
            bed_writer.write(f"{chrom}\t{start}\t{end}\n")
