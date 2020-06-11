# -*- coding: utf-8 -*-

# under command line windows, install required packages.


# How to install python modules without root access?
# https://stackoverflow.com/questions/7465445/how-to-install-python-modules-without-root-access
# pip install --user package_name
# pip install --user pandas
# pip install --user networkx

# If you have access to the root, say, in your own server or laptop.
# python -m pip install pandas
# python -m pip install networkx

import sys
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

import pandas as pd
import networkx as nx

if len(sys.argv) < 3:
    print("Usage: ", sys.argv[0], "inputfile.csv", "outputfile.csv")
    exit(0)

input_file = sys.argv[1]
output_file = sys.argv[2]

# read the edge list to a pandas data frame. The input is the cross weight between nodes.
pd_df = pd.read_csv(input_file)

# Contruct the graph using the pandas data frame edge list of the nodes with cross nodes
nx_graph = nx.from_pandas_edgelist(pd_df, "from", "to", "weight")

# Run the blossom matching algorithm provided by networkx package.

'''
not_maxcardinality_matched = nx.max_weight_matching(nx_graph, maxcardinality=False, weight="weight")
len(not_maxcardinality_matched)

# convert the matching results to a pandas data frame and write to the hard disk.
from_col = [x[0] for x in not_maxcardinality_matched]
to_col = [x[1] for x in not_maxcardinality_matched]
pd_output_df = pd.DataFrame({"from":from_col, "to":to_col}).sort_values("from")
pd_output_df.to_csv(("not_max_cardinality_" + output_file), index=False)
'''

# The same operation except for maxcardinality is set to be true.
maxcardinality_matched = nx.max_weight_matching(nx_graph, maxcardinality=True, weight="weight")
len(maxcardinality_matched)

# convert the matching results to a pandas data frame and write to the hard disk.
from_col = [x[0] for x in maxcardinality_matched]
to_col = [x[1] for x in maxcardinality_matched]
pd_output_df = pd.DataFrame({"from":from_col, "to":to_col}).sort_values("from")
pd_output_df.to_csv(output_file, index=False)
print(input_file, output_file, sep= "\n")
