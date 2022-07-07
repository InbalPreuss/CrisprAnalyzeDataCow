import itertools

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def interleave(iter1, iter2) :
    for (forward, reverse) in itertools.izip(iter1, iter2):
        assert forward.id == reverse.id
        forward.id += "/1"
        reverse.id += "/2"
        yield forward
        yield reverse

    s = SeqIO.SeqRecord(seq='agtc', id='1', name='1', description='1')


def convert_seq_to_seq():

    df = pd.read_excel("C:\\Users\\Inbal\\OneDrive\\Documents\\IDC\\PhD\\CRISPR\\DataVolcaniCow\\btCRISPRko.v1_small_3000.xlsx")

    print(df)
    sequence_set = df['seq']
    sequence_name_set = df['Name']
    records_forward = (SeqRecord(seq=Seq(seq), id=str(name), letter_annotations={'phred_quality': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]}) for seq, name in zip(sequence_set, sequence_name_set))
    SeqIO.write(records_forward, "C:\\Users\\Inbal\\OneDrive\\Documents\\IDC\\PhD\\CRISPR\\DataVolcaniCow\\read1.fq", "fastq")
    records_reverse = (SeqRecord(seq=Seq(seq).reverse_complement(), id=str((index)), letter_annotations={'phred_quality': [1, 2, 3, 4, 5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19, 20]}) for index, seq in enumerate(sequence_set))
    SeqIO.write(records_reverse, "C:\\Users\\Inbal\\OneDrive\\Documents\\IDC\\PhD\\CRISPR\\DataVolcaniCow\\read2.fq", "fastq")


if __name__ == '__main__':
    convert_seq_to_seq()
    # leftReads = SeqIO.parse("/scratch/AiptasiaMiSeq/fastq/Aip02.R1.fastq", "fastq")
    # rightReads = SeqIO.parse("/scratch/AiptasiaMiSeq/fastq/Aip02.R2.fastq", "fastq")

    # handle = open("interleave.fastq", "w")
    # count = SeqIO.write(interleave(leftReads, rightReads), handle, "fastq")
    # handle.close()
    # print("{} records written to interleave.fastq".format(count))