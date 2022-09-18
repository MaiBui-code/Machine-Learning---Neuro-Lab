"""Before running the code, I have to pip install pymaid, navis into my terminal. 

Pymaid (short for “Python-CATMAID”) is a Python library for fetching, 
analyzing and visualizing data generated with CATMAID.

NAVis is a Python library for analysis and visualization of neuron morphology. 
It stands on the shoulders of the excellent natverse for R.

Besides, I also have to create the token as a password to protect the data"""

"""Find the distance of a trachea from the infiltration point to the ends.
Give figures to display and summarize data"""

import navis
import pymaid
import os
import numpy as np 
import pandas as pd
import math 

import matplotlib.pyplot as plt
from fafbseg import flywire

navis.set_pbars(jupyter=False)
pymaid.set_loggers('WARNING')

rm = pymaid.CatmaidInstance("https://spine.itanna.io/catmaid/fafb-edu/#",
                            http_password=None, http_user=None, project_id=4,
                            api_token=os.environ['CATMAID_TOKEN'])         


flywire.set_chunkedgraph_secret(os.environ['FLYWIRE_TOKEN'])

ready = pymaid.get_neuron('annotation:Trachea_Ready to analyze')
ready_list = ready.skeleton_id.tolist()
print(ready_list)

#create a function to have a list of distance from TLC to ends of 
#1 skeleton having id as ids
def distance_tlc_ends(ids):

    #get ends' nodes coordination and turn them into lists
    ends = pymaid.find_nodes(tags="ends", skeleton_ids=ids)
    ends_list_x = ends['x'].tolist()
    ends_list_y = ends['y'].tolist()
    ends_list_z = ends['z'].tolist()

    #get tlc's nodes coordination and turn them into lists
    tlc = pymaid.find_nodes(tags="TLC", skeleton_ids=ids)
    tlc_list_x = tlc['x'].tolist()
    tlc_list_y = tlc['y'].tolist()
    tlc_list_z = tlc['z'].tolist()

    #count the distance from each tlc to ends
    distance = []
    for cor in range(len(ends_list_x)):
        for tlc_cor in range(len(tlc_list_x)):
            x1 = tlc_list_x[tlc_cor]
            y1 = tlc_list_y[tlc_cor]
            z1 = tlc_list_z[tlc_cor]
            x2 = ends_list_x[cor]
            y2 = ends_list_y[cor]
            z2 = ends_list_z[cor]
            dis = math.sqrt((x1 - x2)**2 + (y1 -y2)**2 + (z1 - z2)**2)
            distance.append(dis)

    #return the list of distance 
    return distance

all_distance = []
for ske_id in ready_list:
    distance = distance_tlc_ends(ske_id)
    all_distance.append(distance)

#create a pandas table that have the distance of all the skeletons
# Pandas dataframe
#df = pd.DataFrame(all_distance, index =ready_list,columns =['distance'])
df = pd.DataFrame (all_distance).transpose()
df.columns = ready_list
print (df)

#plot multiple box plots to compare the distance from tlc to ends between skeleton_ids
plt.rcParams["figure.figsize"] = [17.50, 3.50]
plt.rcParams["figure.autolayout"] = True
data = df
ax = data[ready_list].plot(kind='box', title='Distance between TLCs and ends')

# Display the plot
plt.show()
plt.savefig('dis_TLCs_ends')

#create a function to have a list of distance 
#from starters (maybe TLCs/inf/branch) to ends(maybe ends/ uncertain ends/ TLC/...) of 
#1 skeleton having id as ids
def distance_measure(starter, end, ids):

    #get ends' nodes coordination and turn them into lists
    ends = pymaid.find_nodes(tags=end, skeleton_ids=ids)
    ends_list_x = ends['x'].tolist()
    ends_list_y = ends['y'].tolist()
    ends_list_z = ends['z'].tolist()

    #get starter's nodes coordination and turn them into lists
    starters = pymaid.find_nodes(tags=starter, skeleton_ids=ids)
    starters_list_x = starters['x'].tolist()
    starters_list_y = starters['y'].tolist()
    starters_list_z = starters['z'].tolist()

    #count the distance from each tlc to ends
    distance = []
    for cor in range(len(ends_list_x)):
        for starters_cor in range(len(starters_list_x)):
            x1 = starters_list_x[starters_cor]
            y1 = starters_list_y[starters_cor]
            z1 = starters_list_z[starters_cor]
            x2 = ends_list_x[cor]
            y2 = ends_list_y[cor]
            z2 = ends_list_z[cor]
            dis = math.sqrt((x1 - x2)**2 + (y1 -y2)**2 + (z1 - z2)**2)
            distance.append(dis)

    #return the list of distance 
    return distance

    #Loop through all the skeletons to have a list that contains all the distances
all_distance_infil = []
for ske_id in ready_list:
    distance = distance_measure("Infiltration point", "ends", ske_id)
    all_distance_infil.append(distance)

#create a pandas table that have the distance of all the skeletons
# Pandas dataframe
#df = pd.DataFrame(all_distance, index =ready_list,columns =['distance'])
df_infil = pd.DataFrame (all_distance_infil).transpose()
df_infil.columns = ready_list
#plot multiple box plots to compare the distance from tlc to ends between skeleton_ids
plt.rcParams["figure.figsize"] = [12.50, 3.50]
plt.rcParams["figure.autolayout"] = True
data = df_infil
ax = data[ready_list].plot(kind='box', title='Distance between Infiltration points and ends')

# Display the plot
plt.show()
plt.savefig('dis_infil_ends')
print(df_infil)