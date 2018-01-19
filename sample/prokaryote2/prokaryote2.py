import sys 
import numpy as np
sys.path.append("../../")
import pycos as pycos
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
    pycos = pycos.PYCOS()
    pycos.read_locus(record_parse, interspace=0, bottom=300, height=50, plot=False)
    link_dict = read_link(sys.argv[2]) 
    for key,value in link_dict.items():
        pycos.chord_plot(value[0],value[1],top=800,color="#FF7F0E",alpha=0.7) 
    
    locus = list(pycos.locus_dict.keys())[0]
    pycos.calc_cdsdensity(locus,plus=True,minus=False,window_size=10000)
    pycos.calc_cdsdensity(locus,minus=True,plus=False,window_size=10000)
    pycos.plot_data(locus, np.array(pycos.locus_dict[locus]["cds_density_plus"]), bottom=900, height=80, xaxes=False, yaxes=False, plot_style="heatmap") 
    pycos.plot_data(locus, np.array(pycos.locus_dict[locus]["cds_density_minus"]), bottom=800, height=80, xaxes=False, yaxes=False, plot_style="heatmap", cmap=plt.cm.Blues) 
    pycos.save()
