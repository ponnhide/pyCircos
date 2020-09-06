#%matplotlib inline
import sys 
import numpy as np
sys.path.append("../../")
from pycircos import *
from Bio import SeqIO

if __name__ == "__main__":
    gbk = SeqIO.parse("GCF_000005845.2_ASM584v2_genomic.gbff","genbank")
    gcircle = Gcircle()
    gcircle.interspace = 0.0
    gcircle.read_locus(gbk, bottom=900, height=0, linewidth=0) 
    gcircle.set_locus()
    locus_names = list(gcircle.locus_dict.keys())
    gcircle.plot_features('NC_000913.3', bottom=800, height=50, facecolor="#ef6f6a", requirement=lambda x:x.location.strand==-1)
    gcircle.plot_features('NC_000913.3', bottom=750, height=50, facecolor="#6388b4", requirement=lambda x:x.location.strand==1)
    gcircle.save()
