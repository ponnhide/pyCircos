import sys 
import numpy as np
sys.path.append("../../")
from pycircos import * 

if __name__ == "__main__":
    gcircle = Gcircle()
    gcircle.add_locus("1", 1000, bottom=900, height=100) #name, length, bottom (0<=bottom<=1000), height (0<=bottom<=1000)
    gcircle.add_locus("2", 2000, bottom=900, height=100, facecolor="#ED665D")
    gcircle.add_locus("3", 3000, bottom=900, height=100, facecolor="#6DCCDA")
    gcircle.add_locus("4", 2000, bottom=800, height=100, facecolor="#ED97CA")
    gcircle.add_locus("5", 5000, bottom=950, height=50, facecolor="#EDC948")
    gcircle.set_locus() #Creat figure object
  
    data = np.random.rand(100)
    gcircle.scatter_plot("1", data, bottom=900, height=-100)
   
    data = np.random.rand(100)
    gcircle.line_plot("2", data, bottom=600, height=100)
   
    data = np.random.rand(100)
    gcircle.line_plot("2", data, bottom=600, height=100)
    
    data = np.random.rand(50)
    gcircle.heatmap("3", data, bottom=600, height=100)
    
    data = np.random.rand(50)
    gcircle.bar_plot("4", data, bottom=600, height=-200)
    
    data = np.random.rand(50)
    gcircle.bar_plot("4", data, bottom=600, height=200)
 
    gcircle.chord_plot(["4", 1000, 1100], ["1", 200, 400], bottom=400)
    gcircle.chord_plot(["5", 1000, 1100, 950], ["3", 1000, 2000, 600], color="#FF0000")
    gcircle.chord_plot(["5", 4000, 4500, 950], ["2", 500, 1000, 400], color="#F2BE2B")
   
    gcircle.save() 
