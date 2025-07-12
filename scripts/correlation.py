#!/usr/bin/env python3

import pandas as pd
from scipy.stats import pearsonr, spearmanr
import argparse

parser = argparse.ArgumentParser(description='Calculate Pearson or Spearman correlation coefficient between two columns of a dataset.')
parser.add_argument('method', choices=['pearson', 'spearman'], help='Correlation method to use.')
parser.add_argument('dataset', help='Path to the dataset file.')
parser.add_argument('col1', help='Name of the first column to use for correlation.')
parser.add_argument('col2', help='Name of the second column to use for correlation.')

args = parser.parse_args()

data = pd.read_table(args.dataset, sep='\s|\t',engine='python')

column1 = data[args.col1]
column2 = data[args.col2]

if args.method == 'pearson':
    correlation, p_value = pearsonr(column1, column2)
    method_name = 'Pearson'
elif args.method == 'spearman':
    correlation, p_value = spearmanr(column1, column2)
    method_name = 'Spearman'

print(f'{method_name} correlation between {args.col1} & {args.col2}: {correlation}, p-value: {p_value}')

