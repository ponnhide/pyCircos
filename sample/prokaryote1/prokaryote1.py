import sys 
import numpy as np
sys.path.append("../../")
import pycos as pycos
from Bio import SeqIO

if __name__ == "__main__":
    record_parse = SeqIO.parse(sys.argv[1],"genbank")
    pycos = pycos.PYCOS()
    pycos.read_locus(record_parse, interspace=0, bottom=300, height=50, plot=False)
    locus = pycos.locus_dict.keys()[0]
    pycos.calc_gcamount(locus,window_size=1000,slide_size=1000)
    pycos.calc_gcskew(locus,window_size=1000,slide_size=1000)
    pycos.plot_data(locus, np.array(pycos.locus_dict[locus]["gc_amount"]), bottom=350, height=150, xaxes=False, yaxes=False, plot_style="normal", color1="#707070", circular=True)
    pycos.plot_data(locus, np.array(pycos.locus_dict[locus]["gc_skew"]), bottom=580, height=120, xaxes=True, yaxes=False, plot_style="fill")
    pycos.plot_feature(feat_type="gene", bottom=700, requirement=lambda x: x.strand==-1, color="#74AFB9")
    pycos.plot_feature(feat_type="gene", bottom=790, requirement=lambda x: x.strand== 1, color="#C481A1")
    pycos.plot_ticks(bottom=870) 
    pycos.save()
