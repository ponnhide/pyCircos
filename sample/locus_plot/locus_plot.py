import sys 
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
from Bio import SeqIO
sys.path.append("../../")
import pycos as pycos

Set2  = plt.cm.Set2
Set2  = [Set2(i) for i in range(6)]
if __name__ == "__main__":
    record_parse = SeqIO.parse(sys.argv[1],"genbank")
    genome = pycos.GENOME()
    genome.read_locus(record_parse,interspace=0.02, bottom=800, height=80, requirement=lambda x: "NC_0032" in x, color_list=Set2)    
    genome.save(format_type="png")

