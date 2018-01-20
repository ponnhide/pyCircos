import os 
import sys
import math
import collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.path    as mpath
import matplotlib.patches as mpatches
from Bio import SeqIO

###
#This circular plot is for multi locus genbank. YOU can win for circos 
###
tab10 = ["#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","8C564D","#E377C2","#7F7F7F","#BCBD22","#17BECF"]

class GENOME(object): 
    def __init__(self,figsize=(6,6)):
        self.figure = plt.figure(figsize=figsize)
        #Initial Settings, User cannnot touch the following settings 
        self.ax     = plt.subplot(111, polar=True)
        self.ax.set_theta_zero_location("N")
        self.ax.set_theta_direction(-1)
        self.ax.set_ylim(0,1000)
        self.ax.spines['polar'].set_visible(False)
        self.ax.xaxis.set_ticks([])
        self.ax.xaxis.set_ticklabels([])
        self.ax.yaxis.set_ticks([])
        self.ax.yaxis.set_ticklabels([])
        
        self.sum_length = 0
        self.locus_dict  = collections.OrderedDict() 
        self.record_dict = collections.OrderedDict()  
        self.data = [] 
    
    def read_locus(self, record_parse, record_name=None, interspace=0.01, plot=True, bottom=300, height=50, start=0, end=360, 
            lw=1, color_list=["#E3E3E3"], features=True, circular=False, requirement=lambda x: "NC_" in x):
        #The interspace is the space ration between locus  
        if record_name == None:
            record_name = "Record_" + str(len(self.record_dict.keys())) 
        self.record_dict[record_name] = {}
        self.record_dict[record_name]["record"]     = record_parse
        self.record_dict[record_name]["locus_dict"] = {}
        for locus in self.record_dict[record_name]["record"]:
            if  requirement(locus.id) == True:
                self.locus_dict[locus.id] = {} 
                self.locus_dict[locus.id]["length"] = len(str(locus.seq))
                self.locus_dict[locus.id]["seq"]    = str(locus.seq) 
                if features == True: 
                    self.locus_dict[locus.id]["features"] = locus.features 
                self.sum_length += self.locus_dict[locus.id]["length"]
                self.record_dict[record_name]["locus_dict"][locus.id] = self.locus_dict[locus.id]
                
        self.interspace = interspace * (self.sum_length)
        self.sum_length += self.interspace * (len(self.locus_dict.keys()))
        self.theta = np.linspace(0.0, 2 * np.pi, self.sum_length, endpoint=True)
        
        
        s = self.sum_length * (start * 1.0 / 360.0)
        e = 0
        
        if len(color_list) == 1:
            color_list = color_list * len(self.locus_dict.keys()) 
        
        for i, Id in enumerate(self.locus_dict.keys()):
            self.locus_dict[Id]["start"] = int(s) 
            e   = s + (self.locus_dict[Id]["length"] * (end - start) * 1.0 / 360)
            self.locus_dict[Id]["end"]   = int(e)
            pos  = (s+e) * np.pi / self.sum_length
            posl = (s+e) * np.pi / self.sum_length
            post = -(posl + 0.5*np.pi)-np.pi 
            if plot == True:
                self.locus_dict[Id]["bar"] = self.ax.bar([pos,pos], [0,height], bottom=bottom, width=2.0*np.pi*(e-s)/self.sum_length, color=color_list[i], linewidth=1, edgecolor="k")
                self.ax.bar([pos,pos], [0,height], bottom=bottom, width=2.0*np.pi*(e-s)/self.sum_length, color=color_list[i], linewidth=0, edgecolor=color_list[i])
            self.locus_dict[Id]["color"] = color_list[i]
            s = e + self.interspace

    def plot_feature(self, feat_type="repeat_region", bottom=520, height=80, color="#C9BD74", requirement=lambda x: 1):
        #example of requirement is lambda x: "SPADE" in x.qualifiers["note"][0]
        feat_set = set([])
        for Id in self.locus_dict.keys(): 
            locus_info = self.locus_dict[Id]
            for feat in locus_info["features"]:
                feat_set.add(feat.type)
                if feat_type == feat.type and requirement(feat):
                    start = locus_info["start"] + feat.location.start
                    end   = locus_info["start"] + feat.location.end
                    pos   = 2 * start * np.pi / self.sum_length
                    self.ax.bar([pos,pos], [0,height], bottom=bottom, width=2.0*np.pi*(end-start)/self.sum_length, color=color, linewidth=0) 
                
    def calc_gcamount(self, key, window_size=100000, slide_size=100000):
        seq = self.locus_dict[key]["seq"]
        gc_amounts = []
        for i in range(0,len(seq),slide_size):
            gc_amount = (seq[i:i+window_size].upper().count("G") + seq[i:i+window_size].upper().count("C")) * 1.0 / window_size
            gc_amounts.append(gc_amount)
        gc_amounts.append((seq[i:].upper().count("G") + seq[i:i+window_size].upper().count("C")) * 1.0 / (len(seq)-i))
        self.locus_dict[key]["gc_amount"] = gc_amounts
        return gc_amounts
    
    def calc_gcskew(self, key, window_size=100000, slide_size=100000):
        #(G-C)/(G+C) 
        seq = self.locus_dict[key]["seq"]
        gc_skews = []
        for i in range(0,len(seq),slide_size):
            gc_skew = (seq[i:i+window_size].upper().count("G") - seq[i:i+window_size].upper().count("C")) * 1.0 / (seq[i:i+window_size].upper().count("G") + seq[i:i+window_size].upper().count("C")) * 1.0
            gc_skews.append(gc_skew)
        gc_skews.append((seq[i:].upper().count("G") - seq[i:].upper().count("C")) * 1.0 / (seq[i:].upper().count("G") + seq[i:].upper().count("C")) * 1.0)
        self.locus_dict[key]["gc_skew"] = gc_skews
        return gc_skews 

    def calc_cdsdensity(self, key, window_size=100000, plus=True, minus=True):
        gene_num   = 0 
        sum_length = 0
        gene_nums  = [] 
        for feat in self.locus_dict[key]["features"]:
            if feat.type == "CDS" and plus and minus:
                gene_num += 1    
            elif feat.type == "CDS" and feat.strand==1 and plus: 
                gene_num += 1    
            elif feat.type == "CDS" and feat.strand==-1 and minus: 
                gene_num += 1    

            if feat.type == "CDS" and int(feat.location.parts[-1].end) - sum_length > window_size:
                gene_num = gene_num * 1000000 * 1.0/ window_size
                sum_length += window_size
                gene_nums.append(gene_num) 
                gene_num = 0

        gene_num = gene_num * 1000000 * 1.0/(int(feat.location.parts[0].end) - sum_length)
        gene_nums.append(gene_num) 
        if plus and minus:
            self.locus_dict[key]["cds_density"] = gene_nums
        elif plus:
            self.locus_dict[key]["cds_density_plus"] = gene_nums
        else:
            self.locus_dict[key]["cds_density_minus"] = gene_nums
        return gene_nums

    def plot_data(self, key, data, bottom=360, log=False, height=150, xaxes=False, yaxes=False, plot_style="normal", 
            circular=False, lw=0.5, color="k", color1="#D62728", color2="#1F77B4", cmap=plt.cm.Reds):
        #data is composed of # of locus data. It is like [[~],[~]] 
        data = np.array(data)
        if log == True:
            if 0 in data:
                data = np.log10(data+1)
            else:
                data = np.log10(data)
        data = data * height / (np.max(np.abs(data)) - 0) 
        theta = self.theta[self.locus_dict[key]["start"]:self.locus_dict[key]["end"]]

        if len(data) != len(theta):
            new_atheta = [theta[int(len(theta)*1.0/len(data) * j)] for j in range(len(data))]
            data  = data + bottom                
            theta = np.array(new_atheta)
        else:
            data  = data + bottom                
            theta = np.array(atheta)

        if circular == True:
            np.append(data,data[0]) 
            np.append(theta,theta[0]) 

        pos = (self.theta[self.locus_dict[key]["start"]] + self.theta[self.locus_dict[key]["end"]-1]) * 0.5
        if xaxes == True:
            self.ax.bar([pos,pos], [0,3], bottom=bottom, width=self.theta[self.locus_dict[key]["end"]-1] - self.theta[self.locus_dict[key]["start"]], linewidth=0, color="k")
       
        if plot_style == "normal":
            self.ax.plot(theta,data, color="k",lw=lw)
        
        elif plot_style == "fill":  
            self.ax.fill_between(theta,bottom,data,where=data>bottom, facecolor=color1)
            self.ax.fill_between(theta,bottom,data,where=data<bottom, facecolor=color2)

        elif plot_style == "scatter":
            self.ax.sactter(theta,data,where=data>bottom,s=2)
        
        elif plot_style == "heatmap":
            cmaplist = [cmap(i) for i in range(256)]
            width = (self.theta[self.locus_dict[key]["end"]-1] - self.theta[self.locus_dict[key]["start"]])  * 1.0 / data.size
            for i, atheta in enumerate(theta):
                index = 255.0 * (data[i]-min(data)) / (max(data)-min(data)) 
                self.ax.bar([atheta, atheta], [0,height], bottom=bottom, width=width, linewidth=0, color=cmaplist[int(index)]) 
                

    def plot_ticks(self, bottom=900, height=20, width=0.001*np.pi, space=1000000, axes=False):
        for i,key in enumerate(self.locus_dict.keys()):
            locus_info = self.locus_dict[key]
            locus_len  = locus_info["end"]-locus_info["start"] 
            pos = (self.theta[self.locus_dict[key]["start"]] + self.theta[self.locus_dict[key]["end"]-1]) * 0.5
            labelList  = [j for j in range(0,locus_len,space)]
            
            if axes == True:  
                self.ax.bar([pos,pos], [0,2], bottom=bottom, color="k", width=2*np.pi*locus_len/self.sum_length, lw=0) 
            
            for j in range(1,len(labelList)):
                thetal  = 2 * np.pi * (locus_info["start"] + labelList[j]) / self.sum_length
                thetat  = -(thetal + 0.5*np.pi)-np.pi 
                self.ax.bar([thetal,thetal], [0,height], bottom=bottom, color="k", width=width, linewidth=0) #tick
                if thetal > 6:
                    pass 
                #else:
                #    self.ax.text(0.48*np.cos(thetat)+0.5,0.48*np.sin(thetat)+0.5,str(labelList[j]/1000000.0) + "", horizontalalignment='center',rotation=(thetat-0.5*np.pi)*180/np.pi,verticalalignment='center',fontsize=15, transform=self.ax.transAxes) 
        
    def chord_plot(self,start_list, end_list ,top=900, bottom=0, color="#1F77B4", alpha=0.5):
        #srtart_list and end_list is composed of "locus_id", "start", "end". 
        sstart = self.theta[self.locus_dict[start_list[0]]["start"]+start_list[1]] 
        send   = self.theta[self.locus_dict[start_list[0]]["start"]+start_list[2]] 
        
        ostart = self.theta[self.locus_dict[end_list[0]]["start"]+end_list[1]] 
        oend   = self.theta[self.locus_dict[end_list[0]]["start"]+end_list[2]] 
        
        z1 = top - top * math.cos(abs((send-sstart) * 0.5)) 
        z2 = top - top * math.cos(abs((oend-ostart) * 0.5)) 
        if sstart == ostart: 
            pass 
        else:
            Path      = mpath.Path
            path_data = [(Path.MOVETO,  (sstart, top)),
                         (Path.CURVE3,  (sstart, bottom)),     
                         (Path.CURVE3,  (oend,   top)),
                         (Path.CURVE3,  ((ostart+oend)*0.5, top+z2)),
                         (Path.CURVE3,  (ostart, top)),
                         (Path.CURVE3,  (ostart, bottom)),
                         (Path.CURVE3,  (send,   top)),
                         (Path.CURVE3,  ((sstart+send)*0.5, top+z1)),
                         (Path.CURVE3,  (sstart,   top)),
                        ]
            codes, verts = list(zip(*path_data)) 
            path  = mpath.Path(verts, codes)
            patch = mpatches.PathPatch(path, facecolor=color, alpha=alpha, linewidth=0, zorder=0)
            self.ax.add_patch(patch)

    def save(self, file_name="test", format_type="pdf"):
        if format_type == "pdf":
            self.figure.savefig(file_name + ".pdf", bbox_inches="tight")
        else:
            self.figure.savefig(file_name + "." + format_type, bbox_inches="tight", dpi=600)

if __name__ == "__main__":
    record_parse = SeqIO.parse(sys.argv[1],"genbank")
    chuncos = GENOME()
    chuncos.read_locus(record_parse,interspace=0.01,bottom = 850, height=50, requirement=lambda x: "NC_0032" in x)
    chuncos.chord_plot(["NC_003279.8",0,4000000],["NC_003283.11",6000000,10000000],top=500) 
    chuncos.chord_plot(["NC_003280.10",2000000,4000000],["NC_003281.10",4000000,6000000],top=500,color="#FF7F0E") 
    chuncos.chord_plot(["NC_003282.8",2000000,4000000],["NC_003284.9",6000000,8000000],top=500,color="#2CA02C")
    for key in chuncos.locus_dict.keys():
        chuncos.calc_gcamount(key)
        chuncos.calc_gcskew(key)
        chuncos.calc_cdsdensity(key,window_size=50000)
        chuncos.calc_cdsdensity(key,plus=True,minus=False)
        chuncos.calc_cdsdensity(key,minus=True,plus=False)
        chuncos.plot_data(key, np.array(chuncos.locus_dict[key]["cds_density_plus"]), bottom=750, height=50, xaxes=False, yaxes=False, plot_style="heatmap") 
        chuncos.plot_data(key, np.array(chuncos.locus_dict[key]["cds_density_minus"]), bottom=650, height=50, xaxes=False, yaxes=False, plot_style="heatmap", cmap=plt.cm.Blues) 
        chuncos.plot_data(key, np.array(chuncos.locus_dict[key]["cds_density"]), bottom=500, height=100, xaxes=False, yaxes=False, plot_style="fill", color1="#707070")
    chuncos.plot_ticks() 
    chuncos.save()
    #plt.savefig("test.pdf",bbox_inches="tight")
    #plt.savefig(sys.argv[1]+".pdf")


