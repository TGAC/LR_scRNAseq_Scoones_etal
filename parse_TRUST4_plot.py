####
#parse_trust4_output.py
# 2023 06 26
# Script takes the output of TRUST4 that has been run on the MASseq output (whether original or jumpcode versions). This checks 1) the proportion of barcodes possessing chain 1 only, chain 2 only, or both chains 
# and 2) the cell type associated with the barcode. It produces outputs; 1) a list of chain1-only, a list of chain2-only and a list of both chains present, but then also 3 text file lists of cell types; B, abT or gdT barcodes for using with seurat etc. 
# D Wright Earlham Institute 2023
####

#### import necessary modules
import sys
import os
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#### command line inputs 
barcodereport = sys.argv[1] # TRUST4 output file eg sample1A_original_barcode_report.tsv
outstem = sys.argv[2] # naming stem for outputs for any given run of this script 

## outputs to name prior to running for lists... or make in-house from script?
out_Bbarcodes = open(str(outstem+'_Bbarcodes.txt'),'w')
out_abTbarcodes = open(str(outstem+'_abTbarcodes.txt'),'w')
out_gdTbarcodes = open(str(outstem+'_gdTbarcodes.txt'),'w')
out_chain1 = open(str(outstem+'_chain1_only.txt'),'w')
out_chain2 = open(str(outstem+'_chain2_only.txt'),'w')
out_both = open(str(outstem+'_both_chains.txt'),'w')
overview = open(str(outstem+'_parsing_summary.txt'),'w')
#### function declarations

# read in the TRUST4 report and parse appropriately 
def parse_trust4(barcodereport):
        report = open(barcodereport, 'r')
        #outputs
        B_barcodes = [] # list for b-cell barcodes
        abT_barcodes = [] # list for abT cell barcodes
        gdT_barcodes = [] # etc
        chain1_only = []
        chain2_only = []
        both_chains = []

        # For each line in the count file extract details:
        line_cnt = 0 # use a counter for proportions
        for line in report:
            if line.startswith('#'):
                 continue
            array = line.rstrip().split('\t')
            barcode = array[0].strip()
            celltype = array[1].strip()
            chain1 = array[2].strip() # CSV list
            chain2 = array[3].strip() # CSV list
            #second_chain1 = array[4].rstrip() # CSV list
            #second_chain2 = array[5].rstrip() # CSV list
            
            # assign barcode to celltype list
            if celltype == "B":
                 B_barcodes.append(barcode)
            elif celltype == "abT":
                 abT_barcodes.append(barcode)
            elif celltype == "gdT":
                 gdT_barcodes.append(barcode)
            else:
                 print('barcode '+str(barcode)+' has no celltype recorded') # not expected but just in case
            # separate chains into groups and keep a tally, missing data columns have *, careful with the use of * NOT as wildcard!
            if chain1 == "*":
                 chain2_only.append(line)
            if chain2 == "*":
                 chain1_only.append(line)
            if not chain1 == "*" and not chain2 == "*":
                 both_chains.append(line)
            line_cnt = line_cnt+1 # keep track of number of lines in file (number of barcodes/cells)

        # calculate proportion of each type 
        missing1prop = int(len(chain2_only))/int(line_cnt)
        missing2prop = int(len(chain1_only))/int(line_cnt)
        bothprop = int(len(both_chains))/int(line_cnt)
        
        # check proportion of cell type in each chain-file (1, 2, both) and add to report summary  
        chain1_B = 0
        chain1_abT = 0
        chain1_gdT = 0
        chain2_B = 0
        chain2_abT = 0
        chain2_gdT = 0
        both_B = 0
        both_abT = 0
        both_gdT = 0
        for ele in chain1_only:
            array = ele.rstrip().split('\t')
            celltype = array[1].strip()
            if celltype == "B":
                 chain1_B = chain1_B +1
            elif celltype == "abT":
                 chain1_abT = chain1_abT +1
            elif celltype == "gdT":
                 chain1_gdT = chain1_gdT + 1
        for ele in chain2_only:
            array = ele.rstrip().split('\t')
            celltype = array[1].strip()
            if celltype == "B":
                 chain2_B = chain2_B +1
            elif celltype == "abT":
                 chain2_abT = chain2_abT +1
            elif celltype == "gdT":
                 chain2_gdT = chain2_gdT + 1                
        for ele in both_chains:
            array = ele.rstrip().split('\t')
            celltype = array[1].strip()
            if celltype == "B":
                 both_B = both_B +1
            elif celltype == "abT":
                 both_abT = both_abT +1
            elif celltype == "gdT":
                 both_gdT = both_gdT + 1   
        
        # calculate proportions of each total category (both, 1only, 2only) accounted for by each celltype and add to overview
        chain1_only_B = chain1_B / len(chain1_only)
        chain1_only_abT = chain1_abT / len(chain1_only)
        chain1_only_gdT = chain1_gdT / len(chain1_only)
        chain2_only_B = chain2_B / len(chain2_only)
        chain2_only_abT = chain2_abT / len(chain2_only)
        chain2_only_gdT = chain2_gdT / len(chain2_only)
        both_chains_B = both_B / len(both_chains)
        both_chains_abT = both_abT / len(both_chains)
        both_chains_gdT = both_gdT / len(both_chains)

        # calculate proportions of each cell type for each category...this is different from the proportion of each cell type missing each variable as a feature of whole cell count!
        B_bothchains_prop =  both_B / len(B_barcodes)
        abT_bothchains_prop = both_abT / len(abT_barcodes)
        gdT_bothchains_prop = both_gdT / len(gdT_barcodes)
        B_chain1only_prop =  chain1_B / len(B_barcodes)
        abT_chain1only_prop = chain1_abT / len(abT_barcodes)
        gdT_chain1only_prop = chain1_gdT / len(gdT_barcodes)
        B_chain2only_prop =  chain2_B / len(B_barcodes)
        abT_chain2only_prop = chain2_abT / len(abT_barcodes)
        gdT_chain2only_prop = chain2_gdT / len(gdT_barcodes)

        # print out stats from parse to summary file                                
        overview.write("there are "+str(line_cnt)+ " rows in the file "+str(barcodereport)+"\n")
        overview.write("a total of "+str(len(chain2_only))+" ["+str(round(missing1prop,3))+"] are missing chain 1 info, split by: B cells "+"("+str(chain2_B)+" ["+str(round(chain2_only_B, 3))+"]), abT cells "+"("+str(chain2_abT)+" ["+str(round(chain2_only_abT, 3))+"]) and gdT cells "+"("+str(chain2_gdT)+" ["+str(round(chain2_only_gdT, 3))+"])"+"\n")
        overview.write("a total of "+str(len(chain1_only))+" ["+str(round(missing2prop,3))+"] are missing chain 2 info, split by: B cells "+"("+str(chain1_B)+" ["+str(round(chain1_only_B, 3))+"]), abT cells "+"("+str(chain1_abT)+" ["+str(round(chain1_only_abT, 3))+"]) and gdT cells "+"("+str(chain1_gdT)+" ["+str(round(chain1_only_gdT, 3))+"])"+"\n")
        overview.write("a total of "+str(len(both_chains))+" ["+str(round(bothprop,3))+"] have both chains present, split by: B cells "+"("+str(both_B)+" ["+str(round(both_chains_B, 3))+"]), abT cells "+"("+str(both_abT)+" ["+str(round(both_chains_abT, 3))+"]) and gdT cells "+"("+str(both_gdT)+" ["+str(round(both_chains_gdT, 3))+"])"+"\n")
        overview.write("there are "+str(len(B_barcodes))+" B cell barcodes in total in the file"+"\n")
        overview.write("there are "+str(len(abT_barcodes))+" abT cell barcodes in total in the file"+"\n")
        overview.write("there are "+str(len(gdT_barcodes))+" gdT cell barcodes in total in the file")
        
        #write outputs
        #barcode lists
        for ele in B_barcodes:
             out_Bbarcodes.write(str(ele)+'\n')
        for ele in abT_barcodes:
             out_abTbarcodes.write(str(ele)+'\n')
        for ele in gdT_barcodes:
             out_gdTbarcodes.write(str(ele)+'\n')
        # chain-specific rows
        for ele in chain1_only:
             out_chain1.write(str(ele.rstrip())+'\n')
        for ele in chain2_only:
             out_chain2.write(str(ele.rstrip())+'\n')
        for ele in both_chains:
             out_both.write(str(ele.strip())+'\n')
        
        ## PLOTTING
        # plot cell type proportions where each is the proportion of each category for each cell type (so that each cell type has same bar heights in stacked plot)
        bar_plot = pd.DataFrame({'both chains': [B_bothchains_prop, abT_bothchains_prop, gdT_bothchains_prop], 
                                 'chain 1 only': [B_chain1only_prop, abT_chain1only_prop, gdT_chain1only_prop],
                                 'chain 2 only': [B_chain2only_prop, abT_chain2only_prop, gdT_chain2only_prop]},
                                index = ['B ('+str(len(B_barcodes))+')','abT ('+str(len(abT_barcodes))+')','gdT ('+str(len(gdT_barcodes))+')'])
        sns.set(style='ticks', palette='colorblind')
        myplot = bar_plot.plot(kind = 'bar', stacked = True) #, color = ['Blue', 'yellow', 'green'])
        plt.tight_layout()
        plt.xlabel('cell type')
        myplot.set_xticklabels(myplot.get_xticklabels(), rotation = 0) # ensure labels are horizontal
        #myplot.legend(loc='right', bbox_to_anchor=(1, 1)) #move legend to one side
        myplot.legend(loc = 'lower center', bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False)
        plt.ylabel('proportion of cells')
        #plt.title('Ig chain recovery from MASseq data')
        plt.savefig(str(outstem+'_celltype_chain_proportions.pdf'), bbox_inches = 'tight')  

#### Main Program ####
def Main(barcodereport):
     parse_trust4(barcodereport)
Main(barcodereport)
