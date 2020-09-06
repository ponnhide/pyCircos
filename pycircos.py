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

def read_linkdata(f):
    link_dict = collections.defaultdict(list)  
    for line in open(f):
        line   = line.rstrip().split("\t")
        locus  = line[1]
        values = list(map(int,line[2:]))
        link_dict[line[0]].append([locus] + values)
    return link_dict

class Gcircle(object):
    colors = ["#4E79A7","#F2BE2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"]
    cmaps  = [plt.cm.Reds, plt.cm.Blues, plt.cm.Greens, plt.cm.Greys]  
    def __init__(self):
        self.locus_dict = collections.OrderedDict() 
        self.interspace = np.pi / 60 
        self.bottom      = 500 
        self.height      = 50 
        self.facecolor   = "#DDDDDD"
        self.edgecolor   = "#000000"
        self.linewidth   = 1.0 
        self.markersize  = 2.0 
        self.color_cycle = 0 
        self.cmap_cycle  = 0

    def add_locus(self, name, length, bottom=None, height=None, facecolor=None, edgecolor=None, linewidth=None, interspace=None):
        self.locus_dict[name]                 = {}
        self.locus_dict[name]["length"]       = length
        self.locus_dict[name]["features"] = [] 
        
        if bottom is None:
            self.locus_dict[name]["bottom"] = self.bottom
        else:
            self.locus_dict[name]["bottom"] = bottom
        
        if height is None:
            self.locus_dict[name]["height"] = self.height 
        else:
            self.locus_dict[name]["height"] = height

        if facecolor is None: 
            self.locus_dict[name]["facecolor"] = self.facecolor 
        else:
            self.locus_dict[name]["facecolor"] = facecolor

        if edgecolor is None:
            self.locus_dict[name]["edgecolor"] = self.edgecolor
        else:
            self.locus_dict[name]["edgecolor"] = edgecolor

        if interspace is None:
            self.locus_dict[name]["linewidth"] = self.linewidth
        else:
            self.locus_dict[name]["linewidth"] = linewidth

        if interspace is None:
            self.locus_dict[name]["interspace"] = self.interspace
        else:
            self.locus_dict[name]["interspace"] = interspace 


        sum_length       = sum(list(map(lambda x:  self.locus_dict[x]["length"], list(self.locus_dict.keys()))))
        sum_interspace   = sum(list(map(lambda x:  self.locus_dict[x]["interspace"], list(self.locus_dict.keys()))))
        self.theta_list  = np.linspace(0.0, 2 * np.pi - sum_interspace, sum_length, endpoint=True)
        s = 0
        sum_interspace = 0 
        for key in self.locus_dict.keys():
            self.locus_dict[key]["positions"] = sum_interspace + self.theta_list[s:s+self.locus_dict[key]["length"]+1]
            if s+self.locus_dict[key]["length"]+1 > len(self.theta_list):
                self.locus_dict[key]["positions"] = self.locus_dict[key]["positions"] + self.theta_list[:s+self.locus_dict[key]["length"] + 1- len(self.theta_list)]
            s = s + self.locus_dict[key]["length"]
            sum_interspace += self.locus_dict[key]["interspace"]

    def read_locus(self, record_parse, interspace=None, bottom=None,  height=None, facecolor=None, edgecolor=None, linewidth=None, features=True, requirement=lambda x: True):
        #The interspace is the space ration between locus  
        for locus in record_parse:
            if  requirement(locus.id) == True:
                self.locus_dict[locus.id] = {} 
                self.locus_dict[locus.id]["seq"] = str(locus.seq).upper()
                self.locus_dict[locus.id]["length"] = len(str(locus.seq))
                if features == True: 
                    self.locus_dict[locus.id]["features"] = locus.features 
                else:
                    self.locus_dict[locus.id]["features"] = []
        
        if bottom is None:
            bottom = self.bottom
        
        if height is None:
            height = self.height 

        if facecolor is None: 
            facecolor = self.facecolor 

        if edgecolor is None:
            edgecolor = self.edgecolor

        if linewidth is None:
            linewidth = self.linewidth
        
        if interspace is None:
            interspace = self.interspace 

        if type(bottom) != list:
            bottom = [bottom] * len(list(self.locus_dict.keys()))
        
        if type(height) != list:
            height = [height] * len(list(self.locus_dict.keys()))
        
        if type(facecolor) != list:
            facecolor = [facecolor] * len(list(self.locus_dict.keys()))
        
        if type(edgecolor) !=list:
            edgecolor = [edgecolor] * len(list(self.locus_dict.keys()))
        
        if type(linewidth) != list:
            linewidth = [linewidth] * len(list(self.locus_dict.keys()))

        if type(interspace) != list:
            interspace = [interspace] * len(list(self.locus_dict.keys()))
   
        for i, key in enumerate(list(self.locus_dict.keys())):
            self.locus_dict[key]["bottom"]      = bottom[i] 
            self.locus_dict[key]["height"]      = height[i]  
            self.locus_dict[key]["facecolor"]   = facecolor[i]  
            self.locus_dict[key]["edgecolor"]   = edgecolor[i]
            self.locus_dict[key]["linewidth"]   = linewidth[i] 
            self.locus_dict[key]["interspace"]  = interspace[i]  
        
        sum_length       = sum(list(map(lambda x:  self.locus_dict[x]["length"], list(self.locus_dict.keys()))))
        sum_interspace   = sum(list(map(lambda x:  self.locus_dict[x]["interspace"], list(self.locus_dict.keys()))))
        self.theta_list  = np.linspace(0.0, 2 * np.pi - sum_interspace, sum_length, endpoint=True)
        s = 0
        sum_interspace = 0 
        for key in self.locus_dict.keys():
            self.locus_dict[key]["positions"] = sum_interspace + self.theta_list[s:s+self.locus_dict[key]["length"]+1]
            if s+self.locus_dict[key]["length"]+1 > len(self.theta_list):
                self.locus_dict[key]["positions"] = self.locus_dict[key]["positions"] + self.theta_list[:s+self.locus_dict[key]["length"] + 1- len(self.theta_list)]
            s = s + self.locus_dict[key]["length"]
            sum_interspace += self.locus_dict[key]["interspace"]

    def set_locus(self, figsize=(6, 6), lw=1): 
        self.figure = plt.figure(figsize=figsize)
        self.ax     = plt.subplot(111, polar=True)
        self.ax.set_theta_zero_location("N")
        self.ax.set_theta_direction(-1)
        self.ax.set_ylim(0,1000)
        self.ax.spines['polar'].set_visible(False)
        self.ax.xaxis.set_ticks([])
        self.ax.xaxis.set_ticklabels([])
        self.ax.yaxis.set_ticks([])
        self.ax.yaxis.set_ticklabels([])  
                
        pre_e = 0 
        for i, key in enumerate(self.locus_dict.keys()):
            pos       = self.locus_dict[key]["positions"][0] 
            width     = self.locus_dict[key]["positions"][-1] - self.locus_dict[key]["positions"][0]
            height    = self.locus_dict[key]["height"]
            bottom    = self.locus_dict[key]["bottom"]
            facecolor = self.locus_dict[key]["facecolor"]
            edgecolor = self.locus_dict[key]["edgecolor"]
            linewidth = self.locus_dict[key]["linewidth"]
            self.locus_dict[key]["bar"] = self.ax.bar([pos], [height], bottom=bottom, width=width, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor, align="edge")
    
    def add_feature(self, locus_name, start, end, strand, ftype="misc_feature", qualifiers={}):
        if start > len(self.locus_dict[key]["positions"]) - 1:
            raise ValueError("'start' value should be less than length of '{}'".format(locus_name)) 
        
        if end > len(self.locus_dict[key]["positions"]) - 1:
            raise ValueError("'end' value should be less than length of '{}'".format(locus_name)) 

        feat = SeqFeature(FeatureLocation(start, end, strand=strand), type=feature_type)
        for key, value in qualifiers.items():
            if type(value) == list or type(value) == tuple:
                pass
            else:
                qualifiers[key] = [value] 
        
        feat.type = ftype
        feat.qualifiers = qualifiers
        self.locus_dict[locus_name]["features"].append(feat)   
    
    def plot_features(self, key, feat_type="CDS", bottom=500, height=50, facecolor=None, linewidth=0.0, edgecolor="k",  requirement=lambda x: 1):
        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1

        positions = [] 
        widths    = [] 
        locus_info = self.locus_dict[key] 
        for feat in locus_info["features"]:
            if feat_type == feat.type and requirement(feat):
                if feat.location.strand >= 0:
                    s = int(feat.location.parts[0].start.position) 
                    e = int(feat.location.parts[-1].end.position)
                else:
                    s = int(feat.location.parts[-1].start.position) 
                    e = int(feat.location.parts[0].end.position)   

                pos   = locus_info["positions"][s] 
                width = locus_info["positions"][e] - pos    
                positions.append(pos) 
                widths.append(width)
                #self.ax.axvspan(pos, pos+width, bottom, bottom+height, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor)
        
        self.ax.bar(positions, [height] * len(positions), bottom=bottom, width=widths, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor, align="edge") 
       
    def calc_gcratio(self, key, window_size=1000, slide_size=None):
        if slide_size is None:
            slide_size = window_size
        seq = self.locus_dict[key]["seq"]
        gc_amounts = []
        for i in range(0,len(seq),slide_size):
            gc_amount = (seq[i:i+window_size].upper().count("G") + seq[i:i+window_size].upper().count("C")) * 1.0 / window_size
            gc_amounts.append(gc_amount)
        gc_amounts.append((seq[i:].upper().count("G") + seq[i:i+window_size].upper().count("C")) * 1.0 / (len(seq)-i))
        self.locus_dict[key]["gc_ratio"] = gc_amounts
        gc_amounts = np.array(gc_amounts)
        return gc_amounts
    
    def calc_gcskew(self, key, window_size=1000, slide_size=None):
        #(G-C)/(G+C) 
        if slide_size is None:
            slide_size = window_size
        seq = self.locus_dict[key]["seq"]
        gc_skews = []
        for i in range(0,len(seq),slide_size):
            gc_skew = (seq[i:i+window_size].upper().count("G") - seq[i:i+window_size].upper().count("C")) * 1.0 / (seq[i:i+window_size].upper().count("G") + seq[i:i+window_size].upper().count("C")) * 1.0
            gc_skews.append(gc_skew)
        gc_skews.append((seq[i:].upper().count("G") - seq[i:].upper().count("C")) * 1.0 / (seq[i:].upper().count("G") + seq[i:].upper().count("C")) * 1.0)
        self.locus_dict[key]["gc_skew"] = gc_skews
        gc_skews = np.array(gc_skews)
        return gc_skews 

    def calc_feature_density(self, key, feat_type="CDS", window_size=10000, requirement=lambda x: 1):
        sum_length = 0
        value   = 0
        values  = [] 
        for feat in self.locus_dict[key]["features"]:
            if feat.type == feat_type and requirement(feat):
                value += 1
                if int(feat.location.parts[-1].end) - sum_length > window_size:
                    value = value * 1.0/ window_size
                    sum_length += window_size
                    values.append(value) 
                    value = 0 

        value  = value * 1.0/(int(feat.location.parts[0].end) - sum_length)
        values.append(value) 
        return values

    def heatmap(self, locus_name, data, bottom=None, height=None, positions=None, linewidth=0.0, edgecolor="k", cmap=None):
        start = self.locus_dict[locus_name]["positions"][0] 
        end   = self.locus_dict[locus_name]["positions"][-1]
        if positions == None:
            positions = np.linspace(start, end, len(data), endpoint=False)
        else:
            pass 
        
        if height is None:
            height = self.height
        if bottom is None:
            bottom = self.bottom 

        if cmap is None:
            cmap = Gcircle.cmaps[self.cmap_cycle % len(Gcircle.cmaps)] 
            self.cmap_cycle += 1

        width     = positions[1]-positions[0] 
        max_value = max(data) 
        min_value = min(data) 
        bars = self.ax.bar(positions, height=[height] * len(positions), bottom=bottom, align="edge", width=width)  
        for b, bar in enumerate(bars):
            bar.set_facecolor(cmap(data[b]/(max_value-min_value)))
            bar.set_linewidth(linewidth)
            bar.set_edgecolor(edgecolor)
        
    def scatter_plot(self, locus_name, data, bottom=None, height=None, positions=None, markersize=2.0, facecolor=None, linewidth=0.0, edgecolor="k"):
        start = self.locus_dict[locus_name]["positions"][0] 
        end   = self.locus_dict[locus_name]["positions"][-1]
        if positions == None:
            positions = np.linspace(start, end, len(data), endpoint=False)
        else:
            pass 
        
        if bottom is None:
            bottom = self.bottom 

        if height is None:
            top = bottom + self.height
        else:
            top = bottom + height       
        
        if markersize is None:
            markersize = self.markersize
        
        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1

        max_value = max(data) 
        min_value = min(data)
        data = np.array(data) - min_value
        data = bottom + np.array(data * ((top - bottom) / (max_value - min_value)))
        self.ax.scatter(positions, data, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor, s=markersize)
    
    def line_plot(self, locus_name, data, bottom=None, height=None, positions=None, facecolor=None, linewidth=0.5, edgecolor="k", fill=False):
        start = self.locus_dict[locus_name]["positions"][0] 
        end   = self.locus_dict[locus_name]["positions"][-1]
        if positions == None:
            positions = np.linspace(start, end, len(data), endpoint=True)
        else:
            pass 
        
        if bottom is None:
            bottom = seclf.bottom 
        
        if height is None:
            top = bottom + self.height
        else:
            top = bottom + height

        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1

        max_value = max(data) 
        min_value = min(data)
        data = np.array(data) - min_value
        data = bottom + np.array(data * ((top - bottom) / (max_value - min_value)))
        
        if fill==False:
            self.ax.plot(positions, data, color=facecolor, linewidth=linewidth)
        
        elif fill==True:
            if bottom >= 0:
                self.ax.fill_between(positions, data, bottom, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor)
            elif bottom < 0:
                self.ax.fill_between(positions, bottom, data, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor)

    def bar_plot(self, locus_name, data, bottom=None, height=None, positions=None, facecolor=None, linewidth=0.0, edgecolor="k"): 
        start = self.locus_dict[locus_name]["positions"][0] 
        end   = self.locus_dict[locus_name]["positions"][-1]
        if positions == None:
            positions = np.linspace(start, end, len(data), endpoint=False)
        else:
            pass 
        
        if bottom is None:
            bottom = self.bottom 
        
        if height is None:
            top = bottom + self.height
        else:
            top = bottom + height
        
        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1
        
        width = positions[1] - positions[0] 
        max_value = max(data) 
        min_value = min(data)
        data = np.array(data) - min_value
        data = np.array(data * ((top - bottom) / (max_value - min_value)))
        self.ax.bar(positions, data, bottom=bottom, facecolor=facecolor, width=width, linewidth=linewidth, edgecolor=edgecolor, align="edge")

    def tick_plot(self, locus_name, ticks, bottom=800, length=10, spine=False, linewidth=0.5, tickwidth=None, color="k"):
        if tickwidth is None:
            tickwidth = linewidth
        positions = self.locus_dict[locus_name]["positions"]
        for tick in ticks:
            self.ax.plot([positions[tick], positions[tick]], [bottom, bottom+length], color=color, linewidth=tickwidth)
        if spine == True:
            self.ax.plot(positions, [y] * len(positions), color=color, linewidth=linewidth)

    def set_spine(self, locus_name, y=500, linewidth=0.5, color="k"):
        positions = self.locus_dict[locus_name]["positions"]
        self.ax.plot(positions, [y] * len(positions), color=color, linewidth=linewidth)

    def chord_plot(self, start_list, end_list,  bottom=500, center=0, color="#1F77B4", alpha=0.5):
        #start_list and end_list is composed of "locus_id", "start", "end". 
        sstart = self.locus_dict[start_list[0]]["positions"][start_list[1]]
        send   = self.locus_dict[start_list[0]]["positions"][start_list[2]+1]   
        if len(start_list) == 4:
            stop = int(start_list[3]) 
        else:
            stop = bottom

        ostart = self.locus_dict[end_list[0]]["positions"][end_list[1]]
        oend   = self.locus_dict[end_list[0]]["positions"][end_list[2]+1] 
        if len(end_list) == 4:
            etop = int(end_list[3]) 
        else:
            etop = bottom

        z1 = stop - stop * math.cos(abs((send-sstart) * 0.5)) 
        z2 = etop - etop * math.cos(abs((oend-ostart) * 0.5)) 
        if sstart == ostart: 
            pass 
        else:
            Path      = mpath.Path
            path_data = [(Path.MOVETO,  (sstart, stop)),
                         (Path.CURVE3,  (sstart, center)),     
                         (Path.CURVE3,  (oend,   etop)),
                         (Path.CURVE3,  ((ostart+oend)*0.5, etop+z2)),
                         (Path.CURVE3,  (ostart, etop)),
                         (Path.CURVE3,  (ostart, center)),
                         (Path.CURVE3,  (send,   stop)),
                         (Path.CURVE3,  ((sstart+send)*0.5, stop+z1)),
                         (Path.CURVE3,  (sstart, stop)),
                        ]
            codes, verts = list(zip(*path_data)) 
            path  = mpath.Path(verts, codes)
            patch = mpatches.PathPatch(path, facecolor=color, alpha=alpha, linewidth=0, zorder=0)
            self.ax.add_patch(patch)
    
    def save(self, file_name="test", format="pdf", dpi=None):
        self.figure.patch.set_alpha(0.0) 
        if format == "pdf" and dpi is None:
            self.figure.savefig(file_name + ".pdf", bbox_inches="tight")
        else:
            if dpi is None:
                dpi = 600
            self.figure.savefig(file_name + "." + format, bbox_inches="tight", dpi=dpi)
        return self.figure 

if __name__ == "__main__":
    #record_parse = SeqIO.parse(sys.argv[1],"genbank")
    chuncos = Gcircle()
    #chuncos.interspace = np.pi/18
    chuncos.add_locus("1", 1000)
    chuncos.add_locus("2", 2000)
    chuncos.add_locus("3", 3000)
    chuncos.add_locus("4", 2000, bottom=100)
    chuncos.add_locus("5", 5000)
    chuncos.add_locus("5", 5000)

    #chuncos.read_locus(record_parse, requirement=lambda x: "NC_0032" in x, bottom=[400,500,600,700,800,900])
    #chuncos.read_locus(record_parse, requirement=lambda x: "NC_0032" in x)
    
    chuncos.set_locus()
    #chuncos.chord_plot(["NC_003279.8", 0, 4000000, 500], ["NC_003283.11",6000000,10000000, 500]) 
    #chuncos.chord_plot(["NC_003280.10", 2000000, 4000000, 100], ["NC_003281.10", 4000000, 6000000, 500], color="#FF7F0E") 
    #chuncos.chord_plot(["NC_003282.8", 2000000, 4000000], ["NC_003284.9", 6000000, 8000000] ,  color="#2CA02C")

    """
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
    """
    
    chuncos.save()
    #plt.savefig("test.pdf",bbox_inches="tight")


