import math
import numpy as np 
from Bio import Phylo
from pycircos.pycircos import Garc
from pycircos.pycircos import Gcircle 

class Tarc(Garc):
    def __init__(self, arc_id=None, tree=None, format="newick", interspace=3, raxis_range=(900, 950), facecolor=None, edgecolor="#303030", linewidth=0, label=None, labelposition=0, labelsize=10, label_visible=False):
        """
        Parameters
        ----------
        arc_id : str, optional
            Unique identifier for the Garc class object. In the event an id
            value is not provided, an original unique ID is automatically 
            generated for Garc object. The default is None.
        tree : str
            File name of phylogenetic tree
        format : str 
            Format of phylogenetic tree. The default value is "newick."
        interspace : float, optional
            Distance angle (deg) to the adjacent arc section in clockwise 
            sequence. The actual interspace size in the circle is determined by
            the actual arc section width in the resultant circle is determined
            by the ratio of size to the combined sum of the size and interspace
            values of the Garc class objects in the Gcircle class object.
            The default is 3.
        raxis_range : tuple (top=int, bottom=int), optional
            Radial axis range where line plot is drawn.
            The default is (900, 950).
        facecolor : str or tuple representing color code, optional
            Color for filling.. The default is None.
        edgecolor : str or tuple representing color code, optional
            Edge color of the filled area. The default is "#303030".
        linewidth : float, optional
            Edge line width. The default is 0.
        label : str, optional
            Label of the arc section. The default is None.
        labelposition : int, optional
            Relative label height from the center of the arc section.
            The default is 0.
        labelsize : int, optional
            Font size of the label. The default is 10.
        label_visible : bool, optional
            If True, label of the Garc object is shown on the arc section.
            The default is False.

        Returns
        -------
        None
        """
        self.tree = Phylo.read(tree, format)
        self.size = len(self.tree.get_terminals()) 
        self._tree_plotted = 0  
        self._parental_gcircle = None
        if arc_id == None:
            self.arc_id = str(Garc._arcnum) 
        else:
            self.arc_id = arc_id
         
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
        self.labelsize     = labelsize
        
        self._get_col_positions() 
        self._get_row_positions()

        Garc._arcnum += 1
            
    def _get_col_positions(self):
        taxa   = self.tree.get_terminals()
        depths = self.tree.depths()
        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width   = 100 - max_label_width - 1
        if max(depths.values()) == 0:
            depths = tree.depths(unit_branch_lengths=True)
        fudge_margin         = math.log(len(taxa), 2)
        cols_per_branch_unit = (drawing_width - fudge_margin) / float(max(depths.values()))
        positions            = {clade: blen * cols_per_branch_unit + 1.0 for clade, blen in depths.items()}
        self._col_positions   = positions

    def _get_row_positions(self):
        taxa = self.tree.get_terminals()
        positions = {taxon: 2 * idx for idx, taxon in enumerate(taxa)}
        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = (
                positions[clade.clades[0]] + positions[clade.clades[-1]]
            ) // 2
        calc_row(self.tree.root)
        self._row_positions = positions
       
        self.clade_dict       = {} 
        self.terminal_dict    = {} 
        self.nonterminal_dict = {} 
        keys = list(self._row_positions.keys()) 
        keys.sort(key=lambda x: self._row_positions[x]) 
        for key in keys:
            if key in taxa:
                self.terminal_dict[key.name] = key
            else:
                self.nonterminal_dict[key.name] = key
            self.clade_dict[key.name] = key

    def _convert(self, thetalim, rlim):
        col_values = list(self._col_positions.values()) 
        row_values = list(self._row_positions.values()) 
        row_min, row_max = np.min(row_values), np.max(row_values) 
        col_min, col_max = np.min(col_values), np.max(col_values)
        
        theta_positions = [thetalim[0] + (v * abs(thetalim[1]-thetalim[0])/abs(row_min-row_max)) for v in row_values]
        if rlim[0] < rlim[1]:
            r_positions = [rlim[0] + (v * abs(rlim[1]-rlim[0])/abs(col_min-col_max)) for v in col_values]
        else:
            r_positions = [rlim[0] - (v * abs(rlim[1]-rlim[0])/abs(col_min-col_max)) for v in col_values]
        self._theta_dict = dict(zip(list(self._row_positions.keys()) ,theta_positions)) 
        self._r_dict     = dict(zip(list(self._col_positions.keys()) ,r_positions))
    
    def _plot_tree(self, ax, thetalim=None, rlim=None, cladevisual_dict=None, highlight_dict=None, linewidth=None, linecolor=None):
        if linewidth is None:
            linewidth = 0.5
        
        if linecolor is None:
            linecolodr = "k"
        
        if cladevisual_dict is None:
            cladevisual_dict = {} 
        
        self._tree_rlim = rlim
        if rlim[0] < rlim[1]:
            self._tree_direction = "inner"
        else:
            self._tree_direction = "outer"
        
        self._convert(thetalim, rlim)
        s   = []
        c   = []
        ecs = []
        lws = []
        for clade in self._theta_dict:
            if clade.name not in cladevisual_dict:
                cladevisual_dict[clade.name] = {}
            cladevisual_dict[clade.name].setdefault("size",0) 
            cladevisual_dict[clade.name].setdefault("color","k")
            cladevisual_dict[clade.name].setdefault("edgecolor","k")
            cladevisual_dict[clade.name].setdefault("linewidth",0.1)
            s.append(cladevisual_dict[clade.name]["size"])
            c.append(cladevisual_dict[clade.name]["color"])
            ecs.append(cladevisual_dict[clade.name]["edgecolor"])
            lws.append(cladevisual_dict[clade.name]["linewidth"])
        ax.scatter(self._theta_dict.values(), [self._r_dict[clade] for clade in self._theta_dict], s=s, c=c, edgecolors=ecs, linewidths=lws, zorder=1100) 
        for clade in self._theta_dict:
            subclades   = clade.clades
            if len(subclades) > 0:
                sc_thetas   = [self._theta_dict[sc] for sc in subclades] 
                minsc_theta = min(sc_thetas) 
                maxsc_theta = max(sc_thetas) 
                thetas = np.linspace(minsc_theta, maxsc_theta, 100) 
                rs     = [self._r_dict[clade]] * len(thetas) 
                ax.plot(thetas, rs, lw=linewidth, color=linecolor, zorder=0) 
                for sc, sc_theta in zip(subclades, sc_thetas):
                    ax.plot([sc_theta, sc_theta], [self._r_dict[sc], self._r_dict[clade]], lw=linewidth, color=linecolor, zorder=0)
        
        if highlight_dict is not None:
            self.plot_highlight(ax, highlight_dict) 
        self._tree_plotted = 1

    def _plot_highlight(self, ax, highlight_dict=None): 
        if self._tree_plotted == 0:
            raise ValueError("Please run `plot_tree` before running `plot_highlight`") 

        for clade_names in highlight_dict:
            highlight_dict[clade_names].setdefault("color", "#000000")
            highlight_dict[clade_names].setdefault("alpha", 0.25)
            highlight_dict[clade_names].setdefault("label", None)
            highlight_dict[clade_names].setdefault("y", None)
            highlight_dict[clade_names].setdefault("fontsize", self.labelsize)
            
            color    = highlight_dict[clade_names]["color"]
            alpha    = highlight_dict[clade_names]["alpha"]
            label    = highlight_dict[clade_names]["label"]
            fontsize = highlight_dict[clade_names]["fontsize"]
            yloc     = highlight_dict[clade_names]["y"]

            if type(clade_names) is tuple:
                ca     = self.tree.common_ancestor([self.clade_dict[name] for name in clade_names])
                clades = [self.clade_dict[name] for name in clade_names]
            else:
                ca     = self.clade_dict[clade_names]  
                clades = ca.get_terminals() 

            c_thetas   = [self._theta_dict[c] for c in clades] 
            minc_theta = min(c_thetas) 
            maxc_theta = max(c_thetas) 
            
            if self._tree_direction == "inner":
                c_rs       = [self._r_dict[c] for c in clades] 
                maxc_r     = max(c_rs) 
                width = minc_theta - maxc_theta 
                loc   = (minc_theta + maxc_theta) / 2 
                ax.bar([loc], [self._tree_rlim[1] - self._r_dict[ca]], bottom=self._r_dict[ca], width=width, color=color, alpha=alpha, linewidth=0, zorder=1000 + maxc_r)
                if highlight_dict[clade_names]["label"] is None:
                    pass 
                else: 
                    rot = loc*360/(2*np.pi)
                    if 90 < rot < 270:
                        rot = 180-rot
                    else:
                        rot = -1 * rot
                    if yloc is None:
                        yloc = self._tree_rlim[1] - abs(self._tree_rlim[1]-self._tree_rlim[0]) * 0.1
                    ax.text(loc, yloc, str(label), rotation=rot, ha="center", va="center", fontsize=fontsize, zorder=1000 + maxc_r + 0.1)

            else:
                c_rs       = [self._r_dict[c] for c in clades] 
                minc_r     = min(c_rs) 
                width = minc_theta - maxc_theta 
                loc   = (minc_theta + maxc_theta) / 2 
                ax.bar([loc], [self._r_dict[ca] - self._tree_rlim[1]], bottom=self._tree_rlim[1], width=width, color=color, alpha=alpha, linewidth=0, zorder=1000 + (-1 * minc_r))
                if highlight_dict[clade_names]["label"] is None:
                    pass 
                else: 
                    rot = loc*360/(2*np.pi)
                    if 90 < rot < 270:
                        rot = 180-rot
                    else:
                        rot = -1 * rot 
                    if yloc is None:
                        yloc = self._tree_rlim[1] + abs(self._tree_rlim[1]-self._tree_rlim[0]) * 0.1
                    ax.text(loc, yloc, str(label), rotation=rot, ha="center", va="center", fontsize=fontsize, zorder=1000 + (-1 * minc_r) + 0.1)

class Tcircle(Gcircle):
    def __init__(self,  fig=None, figsize=None):
        super().__init__(fig=None, figsize=None)
    
    def __getattr__(self, name):
        if name == "tarc_dict":
            return self._garc_dict

    def add_tarc(self, tarc):
        """
        Add a new Tarc or Garc class object into tarc_dict.

        Parameters
        ----------
        tarc : Tarc class object (default:None)
            Tarc class object to be added.

        Returns
        -------
        None.

        """
        self._garc_dict[tarc.arc_id] = tarc

    def set_tarcs(self, start=0, end=360):
        """
        Visualize the arc rectangles of the Garc class objects in .garc_dict on
        the drawing space. After the execution of this method, a new Garc class
        object cannot be added to garc_dict and figure parameter representing
        maplotlib.pyplot.figure object will be created in Gcircle object.

        Parameters
        ----------
        start : int, optional
            Start angle of the circos plot. The value range is -360 ~ 360.
            The default is 0.
        end : int, optional
            End angle of the circos plot. The value range is -360 ~ 360.
            The default is 360.

        Returns
        -------
        None.
        
        """
        sum_length       = sum(list(map(lambda x:  self._garc_dict[x]["size"], list(self._garc_dict.keys()))))
        sum_interspace   = sum(list(map(lambda x:  self._garc_dict[x]["interspace"], list(self._garc_dict.keys()))))
        start = 2 * np.pi * start / 360
        end   = (2 * np.pi * end / 360) - sum_interspace

        s = 0
        sum_interspace = 0 
        for key in self._garc_dict.keys():
            size = self._garc_dict[key].size
            self._garc_dict[key].coordinates    = [None, None]
            self._garc_dict[key].coordinates[0] = sum_interspace + start + ((end-start) * s/sum_length)
            self._garc_dict[key].coordinates[1] = sum_interspace + start + ((end-start) * (s+size)/sum_length)
            s = s + size
            sum_interspace += self._garc_dict[key].interspace
        
        if self.fig_is_ext:
            self.ax = self.figure.add_axes([0, 0, self.figsize[0], self.figsize[1]], polar=True)
        else:
            self.ax = self.figure.add_axes([0, 0, 1, 1], polar=True)
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
            if facecolor is None:
                facecolor = (0, 0, 0, 0)
            
            if facecolor == (0, 0, 0, 0) and linewidth == 0:
                pass 
            else:
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
    
    def plot_tree(self, tarc_id, rlim=(0,700), cladevisual_dict=None, highlight_dict=None, linecolor="#303030", linewidth=0.5):
        """
        
        Parmeters
        ---------
        tarc_id : str 
            ID of the Tarc class object. The ID should be in Tcircle
            object.tarc_dict.
        rlim : tuple (top=int, bottom=int)
            The top and bottom r limits in data coordinates. 
            The default vlaues is (0, 700).
        cladevisual_dict : dict 
            Dictionay composed of pairs of clade name and a sub-dict holding 
            parameters to visualize the clade. A sub-dict is composed of 
            the following key-value pairs:
                - size: int or float
                    Size of dot. The default value is 5.
                - color: float or str replresenting color code.
                    Face color of dot. The default value is '#303030'.
                - edgecolor: float or str replresenting color code.
                    Edge line color of dot. The default value is '#303030'.
                - linewidth: int or float
                    Edge line width of dot. The default value is 0.5.  
        highlight_dict : dict 
            Dictionay composed of pairs of internal clade name and a sub-dict.
            Instead of clade name, tuples of teminal clade names can also be
            A sub-dict is composed of the following key-value pairs:
                - color: "#000000" 
                    Color of highlight for clades.
                - alpha: flaot 
                    Alpha of highlight for clades. The default vlalue is 0.25.
                - label: str
                    Label. The default vlalue is None.
                - fontsize: int or float 
                    Fontsize of label. The default vlalue is 10.
                - y: int or float
                    Y location of the text. The default vlalue is the bottom 
                    edge of the highliht.
        Returns
        -------
        None
        
        """ 
        
        start      = self._garc_dict[tarc_id].coordinates[0] 
        end        = self._garc_dict[tarc_id].coordinates[-1]
        positions  = np.linspace(start, end, self._garc_dict[tarc_id].size, endpoint=False)
        positions  = positions + abs(positions[1]-positions[0]) * 0.5 
        start, end = positions[0], positions[-1] 
        self._garc_dict[tarc_id]._plot_tree(self.ax, thetalim=(start, end), rlim=rlim, cladevisual_dict=cladevisual_dict, highlight_dict=highlight_dict, linecolor=linecolor, linewidth=linewidth)
    
    def plot_highlight(self, tarc_id, highlight_dict=None):
        """
        tarc_id : str 
            ID of the Tarc class object. The ID should be in Tcircle
            object.tarc_dict.
        highlight_dict : dict 
            Dictionay composed of pairs of internal clade name and a sub-dict.
            Instead of clade name, tuples of teminal clade names can also be
            A sub-dict is composed of the following key-value pairs:
                - color: "#000000" 
                    Color of highlight for clades.
                - alpha: flaot 
                    Alpha of highlight for clades. The default vlalue is 0.25.
                - label: str
                    Label. The default vlalue is None.
                - fontsize: int or float 
                    Fontsize of label. The default vlalue is 10.
                - y: int or float
                    Y location of the text. The default vlalue is the bottom 
                    edge of the highliht.
        
        Returns
        -------
        None

        """
        self._garc_dict[tarc_id]._plot_highlight(self.ax, highlight_dict=highlight_dict)

if __name__ == "__main__":
    pass 
    """
    tree = next(Phylo.parse("../tutorial/sample_data/hmptree.xml", "phyloxml")) 
    col_positions = get_col_positions(tree) 
    row_positions = get_row_positions(tree)
    
    import numpy as np 
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5,5)) 
    ax  = fig.add_axes([0, 0, 1, 1], polar=True)
    theta_dict, r_dict = convert(tree, col_positions, row_positions, (0, np.pi), (0,10))
    
    #ax.scatter(theta_dict.values(), [r_dict[clade] for clade in theta_dict], s=1, color="k") 
    
    plot_lines(ax, theta_dict, r_dict) 

    ax.set_ylim([0,10])
    fig.savefig("test.pdf", bbox_inches="tight")
    """
