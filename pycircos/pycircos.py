import os 
import re
import io 
import sys
import math
import urllib
import tempfile
import requests
import collections
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path    as mpath
import matplotlib.patches as mpatches
from Bio import SeqIO
import Bio

matplotlib.rcParams["figure.max_open_warning"] = 0
matplotlib.rcParams['ps.fonttype']       = 42
matplotlib.rcParams['pdf.fonttype']      = 42
matplotlib.rcParams['font.sans-serif']   = ["Arial","Lucida Sans","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.size']         = 10.0
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0
matplotlib.rcParams['xtick.major.pad']   = 6
matplotlib.rcParams['ytick.major.pad']   = 6
matplotlib.rcParams['xtick.major.size']  = 6
matplotlib.rcParams['ytick.major.size']  = 6

class Garc:
    #list100 = ["#ffcdd2","#f8bbd0","#e1bee7","#d1c4e9","#c5cae9","#bbdefb","#b3e5fc","#b2ebf2","#b2dfdb","#c8e6c9","#dcedc8","#f0f4c3","#fff9c4","#ffecb3","#ffe0b2","#ffccbc","#d7ccc8","#cfd8dc",
    colorlist = ["#ff8a80","#ff80ab","#ea80fc","#b388ff","#8c9eff","#82b1ff","#84ffff","#a7ffeb","#b9f6ca","#ccff90","#f4ff81","#ffff8d","#ffe57f","#ffd180","#ff9e80","#bcaaa4","#eeeeee","#b0bec5",
                 "#ff5252","#ff4081","#e040fb","#7c4dff","#536dfe","#448aff","#18ffff","#64ffda","#69f0ae","#b2ff59","#eeff41","#ffff00","#ffd740","#ffab40","#ff6e40","#a1887f","#e0e0e0","#90a4ae"]
    _arcnum = 0
    def __setitem__(self, key, item):
        self.__dict__[key] = item

    def __getitem__(self, key):
        return self.__dict__[key] 

    def __init__(self, arc_id=None, record=None, size=1000, interspace=3, raxis_range=(500, 550), facecolor=None, edgecolor="#303030", linewidth=0.75, label=None, labelposition=0, labelsize=10, label_visible=False): 
        self._parental_gcircle = None
        if arc_id == None:
            self.arc_id = str(Garc._arcnum) 
        else:
            self.arc_id = arc_id

        if record is None:
            self.record = None
            self.size = size
        
        elif type(record) == Bio.SeqRecord.SeqRecord:
            self.record = record
            self.size   = len(str(self.record.seq))
        
        elif type(record) == str:
            match = re.fullmatch("[a-zA-Z]{1,2}_?[0-9]{5,6}", record)
            if os.path.exists(record) == True:
                self.record = SeqIO.read(value, format="genbank")  
            
            if match is None:
                raise ValueError("Incorrect value for NCBI accession number.") 
            else:
                url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbwithparts&id={}&withparts=on".format(record) 
            outb = io.BytesIO()
            outs = io.StringIO()
            headers = {"User-Agent": "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"} 
            request = urllib.request.Request(url, headers=headers) 
            
            with urllib.request.urlopen(request) as u:
                outb.write(u.read())
            outs.write(outb.getvalue().decode())
            
            with tempfile.TemporaryFile(mode="w+") as o:
                content = outs.getvalue()
                o.write(content)
                o.seek(0)  
                record = SeqIO.parse(o,"genbank")
                record = next(record)
            self.record = record 
            self.size = len(str(self.record.seq))
        else:
            self.record = None
            self.size = size
        
        if facecolor is None:
            facecolor = Garc.colorlist[Garc._arcnum % len(Garc.colorlist)] 
        self.interspace  = 2 * np.pi * (interspace / 360)
        self.raxis_range = raxis_range 
        self.facecolor   = facecolor 
        self.edgecolor   = edgecolor
        self.linewidth   = linewidth
        
        if label is None:
            self.label = arc_id
        else:
            self.label = label

        self.label_visible = label_visible
        self.labelposition = labelposition
        self.labelsize = labelsize
        Garc._arcnum += 1

    def calc_density(self, positions, window_size=1000):
        densities = [] 
        positions.sort()
        for i in range(0, self.size, window_size): 
            source = tuple(range(i, i+window_size))
            amount = 0 
            for pos in positions:
                if type(pos) == int:
                    if pos in source:
                        amount += 1 
                elif type(pos) == tuple:
                    if pos[0] <= source[-1] and pos[1] >= source[0]:
                        amount += 1
                else:
                    raise ValueError("List elements should be int type or tuple consiting of two int values")
            densities.append(amount) 

        source = tuple(range(i,self.size))
        amount = 0 
        for pos in positions:
            if type(pos) == int:
                if pos in source:
                    amount += 1 
            elif type(pos) == tuple:
                if pos[0] <= source[-1] and pos[1] >= source[0]:
                    amount += 1
            else:
                raise ValueError("List elements should be int type or tuple consiting of two int values")
        densities.append(amount*((self.size-i)/window_size))     
        return densities 

    def calc_nnratio(self, n1="G", n2="C", window_size=1000, step_size=None):
        if self.record is None:
            raise ValueError("self.record is None, please specify record value")
        
        if step_size is None:
            step_size = window_size
        
        seq = str(self.record.seq)
        gc_amounts = []
        for i in range(0, len(seq), step_size):
            if n2 is None:
                gc_amount = seq[i:i+window_size].upper().count(n1) * 1.0 / window_size
            else:
                gc_amount = (seq[i:i+window_size].upper().count(n1) + seq[i:i+window_size].upper().count(n2)) * 1.0 / window_size
            gc_amounts.append(gc_amount)
        if n2 is None:
            gc_amounts.append(seq[i:].upper().count(n1) * 1.0 / (len(seq)-i))
        else:
            gc_amounts.append((seq[i:].upper().count(n1) + seq[i:i+window_size].upper().count(n2)) * 1.0 / (len(seq)-i))
        
        self["{}{}_ratio".foramt(n1,n2)] = gc_amounts
        gc_amounts = np.array(gc_amounts)
        return gc_amounts

    def calc_nnskew(self, n1="G", n2="C", window_size=1000, step_size=None):
        #(G-C)/(G+C) 
        if self.record is None:
            raise ValueError("self.record is None, please specify record value")
        
        if step_size is None:
            step_size = window_size
        
        seq = str(self.record.seq) 
        gc_skews = []
        for i in range(0, len(seq), step_size):
            gc_skew = (seq[i:i+window_size].upper().count(n1) - seq[i:i+window_size].upper().count(n2)) * 1.0 / (seq[i:i+window_size].upper().count(n1) + seq[i:i+window_size].upper().count(n2)) * 1.0
            gc_skews.append(gc_skew)
        
        gc_skews.append((seq[i:].upper().count(n1) - seq[i:].upper().count(n2)) * 1.0 / (seq[i:].upper().count(n1) + seq[i:].upper().count(n2)) * 1.0)
        self["{}{}_skew".format(n1,n2)] = gc_skews
        gc_skews = np.array(gc_skews)
        return gc_skews 
    
class Gcircle:
    colors = ["#f44336","#e91e63","#9c27b0","#673ab7","#3f51b5","#2196f3","#00bcd4","#009688","#4caf50","#8bc34a","#cddc39","#ffeb3b","#ffc107","#ff9800","#ff5722","#795548","#9e9e9e","#607d8b"]
    #colors = ["#4E79A7","#F2BE2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"]
    cmaps  = [plt.cm.Reds, plt.cm.Blues, plt.cm.Greens, plt.cm.Greys]  
    
    def __getattr__(self, name):
        if name == "garc_dict":
            return self._garc_dict
    
    def __init__(self, figsize=(8,8)):
        self._garc_dict = {} 
        self.figsize     = figsize
        self.figure      = plt.figure(figsize=figsize)
        self.color_cycle = 0 

    def add_garc(self, garc):
        self._garc_dict[garc.arc_id] = garc

    def set_garcs(self, start=0, end=360):
        sum_length       = sum(list(map(lambda x:  self._garc_dict[x]["size"], list(self._garc_dict.keys()))))
        sum_interspace   = sum(list(map(lambda x:  self._garc_dict[x]["interspace"], list(self._garc_dict.keys()))))
        start = 2 * np.pi * start / 360
        end   = (2 * np.pi * end / 360) - sum_interspace
        #self.theta_list  = np.linspace(0.0, 2 * np.pi - sum_interspace, sum_length, endpoint=True)

        s = 0
        sum_interspace = 0 
        for key in self._garc_dict.keys():
            size = self._garc_dict[key].size
            self._garc_dict[key].coordinates    = [None, None]
            self._garc_dict[key].coordinates[0] = sum_interspace + start + ((end-start) * s/sum_length) #self.theta_list[s:s+self._garc_dict[key]["size"]+1]
            self._garc_dict[key].coordinates[1] = sum_interspace + start + ((end-start) * (s+size)/sum_length)
            #if s+self._garc_dict[key].size+1 > sum_length:
            #    self._garc_dict[key].coordinates[1] = sum_interspace + start + ((end-start) * (s+size+1-sum_length)/sum_length)
            #    self._garc_dict[key].coordinates + self.theta_list[:s+size+1-sum_length]
            s = s + size
            sum_interspace += self._garc_dict[key].interspace
            #print(key, self._garc_dict[key].coordinates) 
        self.figure = plt.figure(figsize=self.figsize)
        self.ax     = plt.subplot(111, polar=True)
        self.ax.set_theta_zero_location("N")
        self.ax.set_theta_direction(-1)
        self.ax.set_ylim(0,1000)
        self.ax.spines['polar'].set_visible(False)
        self.ax.xaxis.set_ticks([])
        self.ax.xaxis.set_ticklabels([])
        self.ax.yaxis.set_ticks([])
        self.ax.yaxis.set_ticklabels([])  
                
        for i, key in enumerate(self._garc_dict.keys()):
            pos       = self._garc_dict[key].coordinates[0] 
            width     = self._garc_dict[key].coordinates[-1] - self._garc_dict[key].coordinates[0]
            height    = abs(self._garc_dict[key].raxis_range[1] - self._garc_dict[key].raxis_range[0])
            bottom    = self._garc_dict[key].raxis_range[0]
            facecolor = self._garc_dict[key].facecolor
            edgecolor = self._garc_dict[key].edgecolor
            linewidth = self._garc_dict[key].linewidth
            #print(key, pos, pos+width) 
            self.ax.bar([pos], [height], bottom=bottom, width=width, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor, align="edge")
            if self._garc_dict[key].label_visible == True:
                rot = (self._garc_dict[key].coordinates[0] + self._garc_dict[key].coordinates[1]) / 2
                rot = rot*360/(2*np.pi)
                if 90 < rot < 270:
                    rot = 180-rot
                else:
                    rot = -1 * rot 
                height = bottom + height/2 + self._garc_dict[key].labelposition
                self.ax.text(pos + width/2, height, self._garc_dict[key].label, rotation=rot, ha="center", va="center", fontsize=self._garc_dict[key].labelsize)
    
    def setspine(self, garc_id, raxis_range=None, facecolor="#30303000", edgecolor="#303030", linewidth=0.75):
        pos     = self._garc_dict[garc_id].coordinates[0] 
        width   = self._garc_dict[garc_id].coordinates[-1] - self._garc_dict[garc_id].coordinates[0]
        height  = abs(raxis_range[1] - raxis_range[0])
        bottom  = raxis_range[0]
        self.ax.bar([pos], [height], bottom=bottom, width=width, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor, align="edge", zorder=0)
    
    def lineplot(self, garc_id, data, positions=None, raxis_range=(550, 600), rlim=None, linestyle="solid", linecolor=None, linewidth=1.0, spine=False):
        start = self._garc_dict[garc_id].coordinates[0] 
        end   = self._garc_dict[garc_id].coordinates[-1]
        size  = self._garc_dict[garc_id].size - 1
        positions_all = np.linspace(start, end, len(data), endpoint=True)
        if positions is None:
            positions = positions_all
        else:
            new_positions = [] 
            for p in positions:
                new_positions.append(start + ((end-start) * p/size))
            positions = new_positions
        
        if raxis_range is None:
            raxis_range = raxis_range[0]   
        bottom = raxis_range[0]
        top    = raxis_range[1] 

        if linecolor is None:
            linecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1
        
        if rlim is None:
            rlim = (min(data) - 0.05 * abs(min(data)), max(data) + 0.05 * abs(max(data))) 

        min_value     = rlim[0]
        max_value     = rlim[1] 
        new_data      = [] 
        new_positions = [] 
        new_data_array      = [] 
        new_positions_array = [] 
        for p, v in zip(positions, data):
            if v > rlim[1] or v < rlim[0]: 
                new_data_array.append(new_data)
                new_positions_array.append(new_positions) 
                new_data      = [] 
                new_positions = [] 
            else: 
                new_data.append(v) 
                new_positions.append(p) 
        new_data_array.append(new_data)
        new_positions_array.append(new_positions) 
        for data, positions in zip(new_data_array, new_positions_array): 
            if len(positions) > 0: 
                data = np.array(data) - min_value
                data = bottom + np.array(data * ((top - bottom) / (max_value - min_value)))
                self.ax.plot(positions, data, color=linecolor, linewidth=linewidth, linestyle=linestyle)
        
        if spine == True:
            self.setspine(garc_id, raxis_range)

    def fillplot(self, garc_id, data, positions=None, raxis_range=(550, 600), rlim=None, base_value=None, facecolor=None, edgecolor="#303030", linewidth=0.0, spine=False):  
        start = self._garc_dict[garc_id].coordinates[0] 
        end   = self._garc_dict[garc_id].coordinates[-1]
        size  = self._garc_dict[garc_id].size - 1
        positions_all = np.linspace(start, end, len(data), endpoint=True)
        if positions is None:
            positions = positions_all
        else:
            new_positions = [] 
            for p in positions:
                new_positions.append(start + ((end-start) * p/size))
            positions = new_positions

        if raxis_range is None:
            raxis_range = raxis_range[0]   
        bottom = raxis_range[0]
        top    = raxis_range[1] 

        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1
        
        if rlim is None:
            rlim = (min(data) - 0.05 * abs(min(data)), max(data) + 0.05 * abs(max(data))) 
        
        min_value     = rlim[0]
        max_value     = rlim[1] 
        if base_value is None:
            base_value = min_value
        new_data      = [] 
        new_positions = [] 
        new_data_array      = [] 
        new_positions_array = [] 
        for p, v in zip(positions, data):
            if v > rlim[1] or v < rlim[0]: 
                new_data_array.append(new_data)
                new_positions_array.append(new_positions) 
                new_data = [] 
                new_positions = [] 
            else: 
                new_data.append(v) 
                new_positions.append(p) 
        new_data_array.append(new_data)
        new_positions_array.append(new_positions) 
        for data, positions in zip(new_data_array, new_positions_array): 
            if len(positions) > 0:  
                base_value = base_value - min_value
                base_value = bottom + base_value * ((top - bottom) / (max_value - min_value)) 
                data = np.array(data) - min_value
                data = bottom + np.array(data * ((top - bottom) / (max_value - min_value)))
                self.ax.fill_between(positions, data, base_value, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor)
        
        if spine == True:
            self.setspine(garc_id, raxis_range)

    def scatterplot(self, garc_id, data, positions=None, raxis_range=(550, 600), rlim=None, markershape="o", markersize=5, facecolor=None, edgecolor="#303030", fadcecolors=None, linewidth=0.0, spine=False):
        start = self._garc_dict[garc_id].coordinates[0] 
        end   = self._garc_dict[garc_id].coordinates[-1]
        size  = self._garc_dict[garc_id].size - 1
        positions_all = np.linspace(start, end, len(data), endpoint=True)
        if positions is None:
            positions = positions_all
        else:
            new_positions = [] 
            for p in positions:
                new_positions.append(start + ((end-start) * p/size))
            positions = new_positions
        
        if raxis_range is None:
            raxis_range = raxis_range[0]   
        bottom = raxis_range[0]
        top    = raxis_range[1] 

        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1
        
        if rlim is None:
            rlim = (min(data) - 0.05 * abs(min(data)), max(data) + 0.05 * abs(max(data))) 

        min_value     = rlim[0]
        max_value     = rlim[1] 
        new_data      = [] 
        new_positions = [] 
        new_data_array      = [] 
        new_positions_array = [] 
        for p, v in zip(positions, data):
            if v > rlim[1] or v < rlim[0]: 
                new_data_array.append(new_data)
                new_positions_array.append(new_positions)
                new_data      = [] 
                new_positions = [] 
            else: 
                new_data.append(v) 
                new_positions.append(p) 
        
        new_data_array.append(new_data)
        new_positions_array.append(new_positions) 
        for positions, data in zip(new_positions_array, new_data_array): 
            if len(positions) > 0:
                data = np.array(data) - min_value
                data = bottom + np.array(data * ((top - bottom) / (max_value - min_value)))
                self.ax.scatter(positions, data, c=facecolor, s=markersize, linewidth=linewidth, edgecolor=edgecolor, marker=markershape)

        if spine == True:
            self.setspine(garc_id, raxis_range)
    
    def barplot(self, garc_id, data, positions=None, width=None, raxis_range=(550, 600), rlim=None, base_value=None, facecolor=None, edgecolor="#303030", linewidth=0.0, spine=False):  
        start = self._garc_dict[garc_id].coordinates[0] 
        end   = self._garc_dict[garc_id].coordinates[-1]
        size  = self._garc_dict[garc_id].size - 1
        positions_all = np.linspace(start, end, len(data), endpoint=True)
        if positions is None:
            positions = positions_all
        else:
            new_positions = [] 
            for p in positions:
                new_positions.append(start + ((end-start) * p/size))
            positions = new_positions
        
        if width is None:
            width = [positions[1] - positions[0]] * len(data) 
        elif type(width) == float or type(width) == int:
            width = [(end-start) * width/size] * len(data)  
        else:
            new_width = [] 
            for w in width:
                new_w = (end-start) * w/size
                new_width.append(new_w) 
            width = new_width 

        if raxis_range is None:
            raxis_range = raxis_range[0]   
        bottom = raxis_range[0]
        top    = raxis_range[1] 

        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1
        
        if rlim is None:
            if min(data) != max(data):
                rlim = (min(data) - 0.05 * abs(min(data)), max(data) + 0.05 * abs(max(data))) 
            else:
                rlim = (min(data), max(data))
        
        min_value     = rlim[0] if rlim[0] is not None else min(data)
        max_value     = rlim[1] if rlim[1] is not None else max(data)
        if base_value is None:
            base_value = min_value

        new_data            = [] 
        new_positions       = [] 
        new_width           = [] 
        new_data_array      = [] 
        new_positions_array = [] 
        new_width_array     = [] 
        for p, v, w in zip(positions, data, width):
            if v > rlim[1] or v < rlim[0]: 
                new_data_array.append(new_data)
                new_positions_array.append(new_positions)
                new_width_array.append(new_width)
                new_data      = [] 
                new_width     = [] 
                new_positions = [] 
            else: 
                new_data.append(v) 
                new_positions.append(p)
                new_width.append(w) 
        
        new_data_array.append(new_data)
        new_positions_array.append(new_positions) 
        new_width_array.append(new_width)
        for data, positions, width in zip(new_data_array, new_positions_array, new_width_array): 
            if len(positions) > 0: 
                base_value = base_value - min_value
                if min_value != max_value:
                    base_value = bottom + base_value * ((top - bottom) / (max_value - min_value)) 
                else:
                    base_value = raxis_range[0] 
                
                data = np.array(data) - min_value
                if min_value != max_value:
                    data = np.array(data) * ((top - bottom) / (max_value - min_value))
                    data = np.array(data) - (base_value - raxis_range[0])
                else:
                    data = [raxis_range[1]-raxis_range[0]] * len(data) 
                self.ax.bar(positions, data, width=width, bottom=base_value, color=facecolor, linewidth=linewidth, edgecolor=edgecolor, align="edge") 
    
        if spine == True:
            self.setspine(garc_id, raxis_range)
    
    def heatmap(self, garc_id, data, positions=None, width=None, raxis_range=(550, 600), cmap=None, vmin=None, vmax=None, edgecolor="#303030", linewidth=0.0, spine=False):  
        start = self._garc_dict[garc_id].coordinates[0] 
        end   = self._garc_dict[garc_id].coordinates[-1]
        size  = self._garc_dict[garc_id].size - 1
        positions_all = np.linspace(start, end, len(data), endpoint=True)
        if positions is None:
            positions = positions_all
        else:
            new_positions = [] 
            for p in positions:
                new_positions.append(start + ((end-start) * p/size))
            positions = new_positions
        
        if width is None:
            width = [positions[1] - positions[0]] * len(data) 
        elif type(width) == float or type(width) == int:
            width = [(end-start) * width/size] * len(data)  
        else:
            new_width = [] 
            for w in width:
                new_w = (end-start) * w/size
                new_width.append(new_w) 
            width = new_width 

        if raxis_range is None:
            raxis_range = raxis_range[0]   
        bottom = raxis_range[0]
        top    = raxis_range[1] 
        height = top - bottom

        if cmap is None:
            cmap = Gcircle.cmaps[self.cmap_cycle % len(Gcircle.cmaps)] 
            self.cmap_cycle += 1

        if vmax is None:
            max_value = max(data)
        else:
            max_value = vmax
        
        if vmin is None:
            min_value = min(data) 
        else:
            min_value = vmin
        
        facecolors = [] 
        for d in data:
            facecolors.append(cmap(d/(max_value-min_value)))
        self.ax.bar(positions, height=[height] * len(positions), width=width, bottom=bottom, color=facecolors, edgecolor=edgecolor, linewidth=linewidth, align="edge")  

        if spine == True:
            self.setspine(garc_id, raxis_range)
    
    def featureplot(self, garc_id, feature_type=None, source=None, raxis_range=(550, 600), facecolor=None, edgecolor="#303030", spine=False):  
        start = self._garc_dict[garc_id].coordinates[0] 
        end   = self._garc_dict[garc_id].coordinates[-1] 
        size  = self._garc_dict[garc_id].size - 1

        if source is None:
            source = self.record.features

        if feature_type is None:
            feature_list = source
        else:
            feature_list = [feat for feat in source if feat.type == feature_type]
        
        positions = [] 
        widths = [] 
        for feat in feature_list:
            if feat.location.strand >= 0:
                s = int(feat.location.parts[0].start.position) 
                e = int(feat.location.parts[-1].end.position)
                pos   = start + ((end-start) * s/size)
                width = start + ((end-start) * e/size) - pos    
                positions.append(pos) 
                widths.append(width)
            else:
                s = int(feat.location.parts[-1].start.position) 
                e = int(feat.location.parts[0].end.position)
                pos   = start + ((end-start) * s/size)
                width = start + ((end-start) * e/size) - pos    
                positions.append(pos) 
                widths.append(width)

        bottom = raxis_range[0]
        top    = raxis_range[1] 
        
        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1
        self.ax.bar(positions, [abs(top-bottom)] * len(positions) , width=widths, bottom=bottom, color=facecolor, align="edge", linewidth=0)
        if spine == True:
            self.setspine(garc_id, raxis_range)
    
    def chord_plot(self, start_list, end_list, facecolor=None, edgecolor=None, linewidth=0.0):
        garc_id1 = start_list[0]
        garc_id2 = end_list[0]
        center = 0 

        start1 = self._garc_dict[garc_id1].coordinates[0] 
        end1   = self._garc_dict[garc_id1].coordinates[-1] 
        size1  = self._garc_dict[garc_id1].size - 1
        sstart = start1 + ((end1-start1) * start_list[1]/size1) 
        send   = start1 + ((end1-start1) * start_list[2]/size1)
        stop   = start_list[3] 
        
        start2 = self._garc_dict[garc_id2].coordinates[0] 
        end2   = self._garc_dict[garc_id2].coordinates[-1] 
        size2  = self._garc_dict[garc_id2].size - 1
        ostart = start2 + ((end2-start2) * end_list[1]/size2) 
        oend   = start2 + ((end2-start2) * end_list[2]/size2)
        etop   = end_list[3] 

        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] + "80" 
            self.color_cycle += 1
        
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
            patch = mpatches.PathPatch(path, facecolor=facecolor, linewidth=linewidth, zorder=0)
            self.ax.add_patch(patch)

    def tickplot(self, garc_id, tick_interval=1000, length=5, width=1.0, direction="outer"):  
        start = self.locus_dict[garc_id].coordinates[0] 
        end   = self.locus_dict[garc_id].coordinates[-1] 
        size  = self.locus_dict[garc_id].size-1
        for i in range(0, size+1, tick_interval):
            for tick in ticks:
                if direction == "outer":
                    self.ax.plot([start + ((end-start) * i/size), start + ((end-start) * i/size)], [bottom, bottom+length], color=color, linewidth=width)  
                if direction == "inner":
                    self.ax.plot([start + ((end-start) * i/size), start + ((end-start) * i/size)], [bottom, bottom-length], color=color, linewidth=width)
   
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
    pass  
