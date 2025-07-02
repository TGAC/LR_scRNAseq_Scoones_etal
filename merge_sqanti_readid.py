import pandas as pd
import optparse

# Parse the command line arguments
parser = optparse.OptionParser()
parser.add_option('-c', '--class_file', dest='file1', help='Path to the first file')
parser.add_option('-s', '--read_id_file', dest='file2', help='Path to the second file')
parser.add_option('-o', '--output_file', dest='output', help='Path to the output file')
parser.add_option('-f', '--filtered_readid', dest='read_id', help='Path to the output filtered read id file')

file1=parser.parse_args()[0].file1
file2=parser.parse_args()[0].file2
outfile=parser.parse_args()[0].output
filtread=parser.parse_args()[0].read_id

df1 = pd.read_csv(file1, sep='\t', usecols=[0, 5, 48], names=['isoform', 'structural_category','filter_result'])
df2 = pd.read_csv(file2, sep='\t', usecols=[0, 1], names=['read_id', 'transcript_id'], compression='gzip')

merged_df = pd.merge(df1, df2, left_on='isoform', right_on='transcript_id')
result_df = merged_df[['isoform', 'structural_category', 'read_id', 'filter_result']]

result_df.to_csv(outfile, sep='\t', index=False)

# filter df by filter_result == 'Isoform' and only keep read_id column
isoform_df = result_df[result_df['filter_result'] == 'Isoform'][['read_id']]
isoform_df.to_csv(filtread, sep='\t', index=False, header=False)
