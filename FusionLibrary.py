import numpy as np

import pandas as pd
from Bio import SeqIO

import sys
sys.path.append('cellar/users/rhahuja/projects/cancer-fusions')
from load_seq_dict import load_CDS_dict, load_refseq_CDS_dict

CDS_dict = load_CDS_dict('/nrnb/users/andreabc/Data/Homo_sapiens.GRCh38.cds.all.fa')
#print(CDS_dict.keys())

fusion_df = pd.read_csv('cosmic/CosmicFusion.tsv', sep="\t")

def reverse(sequence):
    for i in range(len(sequence)):
        new_sequence = sequence[::-1]
        if(new_sequence[i] == "A"):
            new_sequence[i] == "T"
        elif(new_sequence[i] == "C"):
            new_sequence[i] == "G"
        elif(new_sequence[i] == "G"):
            new_sequence[i] == "C"
        elif(new_sequence[i] == "T"):
            new_sequence[i] == "A"
            
        sequence = new_sequence
        
        
required_df = fusion_df.loc[:,["5' Genome start from", "5' Genome stop to","3' Genome start from", "3' Genome stop to", "Translocation Name", "3' Strand", "5' Strand", "Fusion ID"]]
required_df.dropna(inplace = True)
required_df.loc[:,"5' Genome start from"] = required_df["5' Genome start from"].astype(int)
required_df.loc[:,"5' Genome stop to"] = required_df["5' Genome stop to"].astype(int)
required_df.loc[:,"3' Genome start from"] = required_df["3' Genome start from"].astype(int)
required_df.loc[:,"3' Genome stop to"] = required_df["3' Genome stop to"].astype(int)

five_transcript_ID  = required_df.iloc[0,4].split("(")[0].split(".")[0]    
five_start_from     = required_df.iloc[0,0]
five_stop_to        = required_df.iloc[0,1]
five_diff     = five_stop_to - five_start_from
original_sequence_five = CDS_dict[five_transcript_ID]
print(original_sequence_five)

# extract the 1st gene record and find difference between start and end
ofile = open("FusionLibrary.fasta","w")
for i in range(len(required_df)):
    five_transcript_ID  = required_df.iloc[i,4].split("(")[0].split(".")[0]    
    five_start_from     = required_df.iloc[i,0]
    five_stop_to        = required_df.iloc[i,1]
    five_diff     = five_stop_to - five_start_from

#pull sequence for 5' end of new gene
    if(five_transcript_ID in CDS_dict):
        original_sequence_five = CDS_dict[five_transcript_ID]
    else:
        continue
    five_sequence          = original_sequence_five[0:five_diff]
    #print(five_sequence)
    if(required_df.iloc[i,6] == "-"):
        reverse(five_sequence)
        
# extract the 2nd gene record and find difference between start and end
    three_transcript_ID = required_df.iloc[i,4].split("_")
    #print(three_transcript_ID)
    if(len(three_transcript_ID) < 3):
        continue
    else:
        three_transcript_ID = three_transcript_ID[2].split("(")[0].split(".")[0]
    three_start_from    = required_df.iloc[i,2]
    three_stop_to       = required_df.iloc[i,3]
    three_diff    = three_stop_to - three_start_from

#pull sequence for 3' end of new gene
    if(three_transcript_ID in CDS_dict):
        original_sequence_three = CDS_dict[three_transcript_ID]
    else:
        continue
    three_sequence          = original_sequence_three[-1:-1 - three_diff] 
    if(required_df.iloc[i,5] == "-"):
        reverse(three_sequence)
#new sequence created
    new_sequence = five_sequence + three_sequence
   # print(new_sequence)

#generating the header for the sequence in the fasta file

    header = ">" + str(required_df.iloc[i,7]) + "|" + "Gene ENSTs: " + required_df.iloc[i,4].split("(")[0] + required_df.iloc[i,4].split("_")[2].split("(")[0] 
    #print(header)

#generating the fasta file

    ofile.write(header + "\n" + str(new_sequence) + "\n")

ofile.close()



