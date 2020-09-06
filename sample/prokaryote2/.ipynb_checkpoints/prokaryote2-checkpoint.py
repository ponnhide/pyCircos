import sys 
import numpy as np
sys.path.append("../../")
import pycircos as pycircos
from Bio import SeqIO
import collections 
import matplotlib.pyplot as plt

def read_link(f):
    link_dict = collections.defaultdict(list)  
    for line in open(f):
        line = line.rstrip().split("\t")
        link_dict[line[0]].append([line[1],int(line[2]),int(line[3])])
    return link_dict

if __name__ == "__main__":
    record_parse = SeqIO.parse(sys.argv[1],"genbank")
    genome = pycircos.GENOME()
    genome.read_locus(record_parse, interspace=0, bottom=300, height=50, plot=False)
    link_dict = read_link(sys.argv[2]) 
    for key,value in link_dict.items():
        genome.chord_plot(value[0],value[1],top=800,color="#FF7F0E",alpha=0.7) 
    
    locus = list(genome.locus_dict.keys())[0]
    genome.calc_cdsdensity(locus,plus=True,minus=False,window_size=5000)
    genome.calc_cdsdensity(locus,minus=True,plus=False,window_size=5000)
    genome.plot_data(locus, np.array(genome.locus_dict[locus]["cds_density_plus"]), bottom=900, height=100, xaxes=False, yaxes=False, plot_style="heatmap") 
    genome.plot_data(locus, np.array(genome.locus_dict[locus]["cds_density_minus"]), bottom=800, height=100, xaxes=False, yaxes=False, plot_style="heatmap", cmap=plt.cm.Blues) 
    genome.save()
