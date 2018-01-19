import sys 
import collections
import numpy as np
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
sys.path.append("../../")
import pycos as pycos
from Bio import SeqIO

#tab6 = ["#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564D"] 
Set2  = plt.cm.Set2
Set2  = [Set2(i) for i in range(6)]

def read_link(f):
    link_dict = collections.defaultdict(list)  
    for line in open(f):
        line = line.rstrip().split("\t")
        link_dict[line[0]].append([line[1],int(line[2]),int(line[3])])
    return link_dict


if __name__ == "__main__":
    record_parse = SeqIO.parse(sys.argv[1],"genbank")
    pycos = pycos.PYCOS()
    pycos.read_locus(record_parse,interspace=0.02, bottom=790, height=30, requirement=lambda x: "NC_0032" in x, color_list=Set2)
    
    link_dict  = read_link(sys.argv[2]) 
    locus_keys = pycos.locus_dict.keys()  
    for key,value in link_dict.items():
        pycos.chord_plot(value[0],value[1],top=700, color=Set2[locus_keys.index(value[0][0])] ,alpha=0.8) 
    
    for key in pycos.locus_dict.keys():
        pycos.calc_gcamount(key)
        pycos.calc_gcskew(key)
        pycos.calc_cdsdensity(key,window_size=50000)
        pycos.calc_cdsdensity(key,plus=True,minus=False)
        pycos.calc_cdsdensity(key,minus=True,plus=False)
        pycos.plot_data(key, np.array(pycos.locus_dict[key]["cds_density_plus"]), bottom=650,  height=60, xaxes=False, yaxes=False, plot_style="heatmap") 
        pycos.plot_data(key, np.array(pycos.locus_dict[key]["cds_density_minus"]), bottom=720, height=60, xaxes=False, yaxes=False, plot_style="heatmap", cmap=plt.cm.Blues) 
        pycos.plot_data(key, np.array(pycos.locus_dict[key]["cds_density"]), bottom=830, height=100, xaxes=False, yaxes=False, plot_style="fill", color1="#777777") 

    pycos.plot_ticks(axes=True, bottom=950) 
    pycos.save()

