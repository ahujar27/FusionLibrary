#from collections import defaultdict
from Bio import SeqIO

# load CDS/cDNA file
def load_CDS_dict(CDS_path='/nrnb/users/andreabc/Data/Homo_sapiens.GRCh38.cds.all.fa'):
    CDS_dict = dict()
    with open(CDS_path, 'r') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            enst = seq_record.id.split('.')[0]
            if enst in CDS_dict:
                # check which one is more accurate
                print(enst) #### none are duplicates!
            if str(seq_record.seq) == 'Sequenceunavailable':
                continue
            CDS_dict[enst] = str(seq_record.seq)  
    return CDS_dict

def load_refseq_CDS_dict(CDS_path='/nrnb/users/andreabc/Data/GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna'):
    CDS_dict = defaultdict(str)
    with open(CDS_path, 'r') as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            np_str = record.id.split('_cds_')[1]
            CDS_dict[np_str] = str(record.seq)
    return CDS_dict

def load_mouse_CDS_dict(CDS_path='/nrnb/users/andreabc/Data/Mus_musculus.GRCm38.cds.all.fa'):
    CDS_dict = defaultdict(list)
    with open(CDS_path, 'r') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            enst = seq_record.id.split('.')[0]
            if enst in CDS_dict:
                # check which one is more accurate
                print(enst) #### none are duplicates!
            if str(seq_record.seq) == 'Sequenceunavailable':
                continue
            CDS_dict[enst] = str(seq_record.seq)  
    return CDS_dict


def load_peptide_dict(path='/nrnb/users/andreabc/Data/Homo_sapiens.GRCh38.pep.all.fa'):
    peptide_dict = defaultdict(list)
    with open(path, 'r') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            enst = seq_record.description.split('transcript:')[1].split('.')[0]
            if enst in peptide_dict:
                # check which one is more accurate
                print(enst) #### none are duplicates!
            if str(seq_record.seq) == 'Sequenceunavailable':
                continue
            peptide_dict[enst] = str(seq_record.seq)  
    return peptide_dict

