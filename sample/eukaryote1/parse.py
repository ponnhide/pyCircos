import os
import sys

n = 0 
flag = 0

for line in open(sys.argv[1]):
    if line[0] == "a":
        if int(line.rstrip().split(" ")[1].split("=")[-1]) >= 1000:
            flag = 1 
            n += 1   
    
    elif line[0] == "s" and flag == 1:
        line = line.replace("  "," ") 
        elements = line.rstrip().split(" ")
        start    = int(elements[-5]) 
        end      = start + int(elements[-4])
        flag    += 1
        print "segdup_" + str(n), elements[1], start, end 
    
    elif line[0] == "s" and flag == 2:
        line = line.replace("  "," ") 
        elements = line.rstrip().split(" ")
        start    = int(elements[-5]) 
        end      = start + int(elements[-4])
        print "segdup_" + str(n), elements[1], start, end 

    else:
        flag = 0

