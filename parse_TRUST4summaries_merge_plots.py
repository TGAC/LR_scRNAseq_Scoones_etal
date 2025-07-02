####
# parse_trust4summaries_merge_plots.py
# 2023 12 06
# Script takes list of parsing summaries (produced by parse_TRUST4_plot.py script) and filenaming stem. It produces an IgG chains proportion plot of all the TRUST4 runs merged
# D Wright / chatGPT used to aid Altair encoding
####

#### import necessary modules
#from calendar import c
import sys
#import seaborn as sns
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import altair as alt

#### command line inputs for script to function
summary_files = sys.argv[1] # file list (full paths) of summary_file.txt for each of the TRUST4 parsing outputs 
outstem = sys.argv[2] # naming stem for outputs for any given run of this script 

def flatten_dict(d, parent_key='', sep='_'):
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def prep_df(df, name):
    df = df.stack().reset_index()
    df.columns = ['runID', 'chains', 'proportion']
    df['celltype'] = name
    return df

def create_merged_barplot(summary_files):
    # dictionary for all the summary data
    # print(summaries)
    prop_dict = {} 
    with open(summary_files) as filelist:
        for file in filelist:
            filename = file.rstrip()
            #print(filename)
            runID=filename.split('/')[-1].split('_defaultsettings')[0].strip() # extract the experiment from the filepath naming convention
            #print(runID)
            prop_dict[runID] = {} # assigning experiment as initial level key
            summary = open(filename, 'r')
            for line in summary: # pull out the relevant info for plotting from the summar lines (all files identical)
                if 'B cell barcodes in total' in line: 
                    Bcells = line.split('there are ')[-1].split(' B cell')[0].strip()
                elif 'abT cell barcodes in total' in line: 
                    abTcells = line.split('there are ')[-1].split(' abT cell')[0].strip()
                elif 'gdT cell barcodes in total' in line:
                    gdTcells = line.split('there are ')[-1].split(' gdT cell')[0].strip() 
                elif 'missing chain 1' in line:
                    Bchain2only = line.split('B cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                    abTchain2only = line.split('abT cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                    gdTchain2only = line.split('gdT cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                elif 'missing chain 2' in line:
                    Bchain1only = line.split('B cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                    abTchain1only = line.split('abT cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                    gdTchain1only = line.split('gdT cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                elif 'both chains present':
                    Bbothchains = line.split('B cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                    abTbothchains = line.split('abT cells ')[-1].split(' [')[0].split('(')[-1].strip() 
                    gdTbothchains = line.split('gdT cells ')[-1].split(' [')[0].split('(')[-1].strip() 
            #print(Bcells)
            #print(Bchain2only)
            #print(Bchain1only)
            #print(Bbothchains)
            
            # calculate plotting proportions as per the other script (simply as a checkpoint - should be the same as previous script output!)
            B_bothchains_prop =  float(Bbothchains) / float(Bcells)
            abT_bothchains_prop =  float(abTbothchains) / float(abTcells)
            gdT_bothchains_prop = float(gdTbothchains)  / float(gdTcells)
            B_chain1only_prop =  float(Bchain1only)  / float(Bcells)
            abT_chain1only_prop = float(abTchain1only) / float(abTcells)
            gdT_chain1only_prop = float(gdTchain1only) / float(gdTcells)
            B_chain2only_prop =  float(Bchain2only) / float(Bcells)
            abT_chain2only_prop = float(abTchain2only) / float(abTcells)
            gdT_chain2only_prop = float(gdTchain2only) / float(gdTcells)
            #print(str(B_bothchains_prop))
            #print(str(B_chain1only_prop))
            #print(str(B_chain2only_prop))
            prop_dict[runID]['B-cells'] = {'Bothchains': B_bothchains_prop, 'chain1only': B_chain1only_prop, 'chain2only': B_chain2only_prop}
            prop_dict[runID]['abT-cells'] = {'Bothchains': abT_bothchains_prop, 'chain1only': abT_chain1only_prop, 'chain2only': abT_chain2only_prop} 
            prop_dict[runID]['gdT-cells'] = {'Bothchains': gdT_bothchains_prop, 'chain1only': gdT_chain1only_prop, 'chain2only': gdT_chain2only_prop} 
            #print(prop_dict)
    
    # Flatten the nested dictionary
    flattened_dict = {k: flatten_dict(v) for k, v in prop_dict.items()}

    # Convert flattened dictionary to a DataFrame
    mergeddata = pd.DataFrame.from_dict(flattened_dict, orient='index')
    mergeddata = mergeddata.rename_axis('experimentID')
    #print(mergeddata)

    # split dataframe to prep/melt for stacked grouping with Altair
    df_B = mergeddata.iloc[:, 0:3].copy()
    df_abT = mergeddata.iloc[:, 3:6].copy()    
    df_gdT = mergeddata.iloc[:, 6:9].copy()
    #print(df_B)
    #print(df_abT)
    #print(df_gdT)
    df_B = prep_df(df_B, 'B')
    df_abT = prep_df(df_abT, 'abT')
    df_gdT = prep_df(df_gdT, 'gdT')
    df = pd.concat([df_B, df_abT, df_gdT])
    print(df)
    #rephrase the label fields for manuscript consistency and clarity
    new_index_labels = ['ONT','ONT','ONT','Revio+JC','Revio+JC','Revio+JC','Revio','Revio','Revio','Sequel-iie+JC','Sequel-iie+JC','Sequel-iie+JC','Sequel-iie','Sequel-iie','Sequel-iie']*3
    new_chain_labels = ['both chains','chain 1 only','chain 2 only']*15
    df['runID'] = new_index_labels
    df['chains'] = new_chain_labels
    print(df)    
    
    #### PLOTTING ####
    #FYI seaborn colorblind palette hex codes: ['#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9']  
    # Plot grouped stacked barplot of all 5 TRUST4 runs with ALTAIR
    figure = alt.Chart(df).mark_bar().encode(
        # tell Altair which field to group columns on
        x=alt.X('celltype:N', axis=alt.Axis(title=None)),
        # tell Altair which field to use as Y values and how to calculate
        y=alt.Y('sum(proportion):Q',
            axis=alt.Axis(
                grid=False,
                title='Proportion of Cells')),
        # tell Altair which field to use to use as the set of columns to be  represented in each group
        column=alt.Column('runID:N', title=None),
        # tell Altair to reorder the columns to have blue on the both-chains (blue bar) at base of stacks
        order=alt.Order('chains:N', sort='ascending'), 
        # tell Altair which field to use for color segmentation 
        color=alt.Color('chains:N', title='IgG Chains',
                scale=alt.Scale(
                    # Chose to use blue, yellow and green as per original plots
                    range=['#0173b2', '#de8f05','#029e73'],
                ),
            ))\
        .configure_view(
            # remove grid lines around column clusters
            strokeOpacity=0
        )\
        .configure_axisX(labelAngle = 0).configure_axis(labelPadding=0.2)
    figure.save(str(outstem+'_celltype_chain_proportions_merged.pdf'))
    figure.save(str(outstem+'_celltype_chain_proportions_merged.png'))

#### Main Program ####
def Main(summary_files):
    create_merged_barplot(summary_files)
Main(summary_files)
