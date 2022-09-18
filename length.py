"""Before running the code, I have to pip install pymaid, navis into my terminal. 

Pymaid (short for “Python-CATMAID”) is a Python library for fetching, 
analyzing and visualizing data generated with CATMAID.

NAVis is a Python library for analysis and visualization of neuron morphology. 
It stands on the shoulders of the excellent natverse for R.

Besides, I also have to create the token as a password to protect the data"""

"""Find the length of a trachea from the infiltration point to the nearest TLC.
And compare them with the total length of the trachea. 
Give figures to display and summarize data"""

import navis
import pymaid
import os
import numpy as np 
import pandas as pd
import math

import matplotlib.pyplot as plt
from fafbseg import flywire
#fafbseg is a set of Python tools to work with various kinds of segmentation data
#in the FAFB dataset: Google's auto-segmentation.

navis.set_pbars(jupyter=False)
pymaid.set_loggers('WARNING')

rm = pymaid.CatmaidInstance("https://spine.itanna.io/catmaid/fafb-edu/#",
                            http_password=None, http_user=None, project_id=4,
                            api_token=os.environ['CATMAID_TOKEN'])         


flywire.set_chunkedgraph_secret(os.environ['FLYWIRE_TOKEN'])

#create a list of trachea ID ready to analyze called ready_list
ready = pymaid.get_neuron('annotation:Trachea_Ready to analyze')
ready_list = ready.skeleton_id.tolist()

#create a list of the total length of a trachea
cable_length = ready.cable_length.tolist()
#create a node list that conatins all the nodes (with specified tag) of a trachea 
def find_all_nodes_list(sk_id, tag):
    node = pymaid.find_nodes(tags=tag, skeleton_ids=sk_id)
    node_list = node.node_id.tolist()
    return node_list
def find_len(start, end, table_):
    #create a function to find the parent id of a given id
    def find_parent(id):
        index = table_[table_['node_id']== id].index.values[0]
        parent = table_.loc[index][1]
        return parent
    
    #find the parent id for a given node and add it to a list
    parent_list = []
    id = end
    while id != start:
        parent_id = find_parent(id)
        parent_list.append(parent_id)
        id = parent_id
    #find the coordinates from the nodes 
    filtered_table = table_.loc[table_['node_id'].isin(parent_list)]
    coor_table = filtered_table[['node_id','x','y','z']]
    
    #rearrange the order of the table so that it follows the order of the parent list
    index_list = []
    for id_index in parent_list:
        index = coor_table[coor_table['node_id']== id_index].index.values[0]
        index_list.append(index)
    reorder_table = coor_table.reindex(index_list)
    
    #change all of the coordinations into 3 lists: x, y, z
    list_x = reorder_table['x'].tolist()
    list_y = reorder_table['y'].tolist()
    list_z = reorder_table['z'].tolist()

    #count the distance (goal)
    total_dis = 0
    for i in range(len(list_x) - 1):
        dis = math.sqrt((list_x[i] - list_x[i+1])**2 + (list_y[i] -list_y[i+1])**2 + (list_z[i] - list_z[i+1])**2)
        total_dis += dis 
    return total_dis
    return start_list, end_list
def find_all_length(sk_id, tag1, tag2):
    list1 = find_all_nodes_list(sk_id, tag1)
    list2 = find_all_nodes_list(sk_id, tag2)
    length_list = []
    end_list = []
    #create a node table of a skeleton id
    table_ = pymaid.get_node_table(sk_id)
    for start_id in list1:
        for end_id in list2:
            try:
                length = find_len(start_id,end_id, table_)
                length_list.append(length)
                end_list.append(end_id)
            except:
                pass
    return length_list

#try the function on the first trachea
len_answer= find_all_length(3244129,"Infiltration point","TLC")
print(len_answer)
print(min(len_answer))

#try another trachea
len_answer= find_all_length(5592325,"Infiltration point","TLC")
print(len_answer)
print(min(len_answer))

all_len = []
for sk_id in ready_list:
    len_list = find_all_length(sk_id,"Infiltration point","TLC")
    all_len.append(len_list)
df_len_tlc = pd.DataFrame (all_len).transpose()
print(df_len_tlc)

df_len_tlc.columns = ready_list
#plot multiple box plots to compare the distance from tlc to ends between skeleton_ids
plt.rcParams["figure.figsize"] = [18.50, 3.50]
plt.rcParams["figure.autolayout"] = True
data = df_len_tlc
ax = data[ready_list].plot(kind='box', title='Lengthe between Infiltration points and TLCs')

# Display the plot
plt.show()
plt.savefig('len_infil_tlc')
print(df_len_tlc)

#find all length from infil to nearest tlc
all_nearest = []
for index,len_list in enumerate(all_len):
    #if the list is not empty (means if the trachea has tlc), 
    #count the min and add it to the list
    if len(len_list) > 0:
        nearest = min(len_list)
        all_nearest.append(nearest)
    else:
        print(index)
        all_nearest.append(-100)
all_nearest

mean = sum(all_nearest)/len(all_nearest)

# Pandas dataframe
#df = pd.DataFrame(all_distance, index =ready_list,columns =['distance'])
df_nearest = pd.DataFrame (all_nearest).transpose()
df_nearest.columns = ready_list
#plot multiple box plots to compare the distance from tlc to ends between skeleton_ids
plt.rcParams["figure.figsize"] = [18.50, 3.50]
plt.rcParams["figure.autolayout"] = True
data = df_nearest
ax = data[ready_list].plot(kind='box', title='Distance between Infiltration points and ends')

# Display the plot
plt.show()
plt.savefig('len_min_tlc')
print(df_nearest)

ratio = []
for index in range(len(all_nearest)):
    tr_ratio = all_nearest[index]/cable_length[index]
    ratio.append(tr_ratio)
fig = plt.figure(figsize =(16, 7))
plt.boxplot(ratio)
plt.show
