#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import csv
threshold = 0.5

output_filename = "output_cluster_based.csv"
data_file = "20190816.Goldstein.integrated5_2_min3.filt.RNA.counts.txt"
cell_cluster_correlation_file = "20190816.Goldstein.integrated5_2_min3.meta.data.txt"

def write_to_file(row):
    with open(output_filename, mode='a') as file:
        file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(row)


# In[4]:


### MAKE LIST OF ALL CELLS IN NEURON CLUSTER     cells_in_neuron_cluster = [[name1, cluster1], [...]


cells_in_neuron_cluster = []
with open(cell_cluster_correlation_file) as f:
    for x, i in enumerate(f):
        line_split = i.split("\n")[0].split("\t")
        if line_split[-2] in ["34", "48"]:
            try:
                cells_in_neuron_cluster.append([line_split[0], int(line_split[-2])]) #CellID, cluster
            except: 
                print(line_split)
### make index list of cells in datafile
cell_indices_in_datafile = []

with open(data_file) as f:
    for x, i in enumerate(f):
        line_split = i.split("\n")[0].split("\t")
        if x == 0:
            list_of_cells = line_split
            break

for i in cells_in_neuron_cluster:
    cell_indices_in_datafile.append(list_of_cells.index(i[0]))


# In[3]:


### ADD ORs for each cell           cells_in_neuron_cluster = [[name1, cluster1, nr_ORs, ORnames]

if os.path.exists(output_filename):
    os.remove(output_filename)

for x, i in enumerate(cells_in_neuron_cluster):
    cells_in_neuron_cluster[x].append(0)

with open(data_file) as f:
    for x, i in enumerate(f):
        line_split = i.split("\n")[0].split("\t")
        if line_split[0][:2] == "OR" or line_split[0][:5] == "VN1R1":
            if line_split[0][2].isdigit() == True:
                for p, c in enumerate(cell_indices_in_datafile):
                    if float(line_split[c+1]) >= threshold:
                        cells_in_neuron_cluster[p][2] += 1
                        cells_in_neuron_cluster[p].append(line_split[0])
                        cells_in_neuron_cluster[p].append(line_split[c+1])
        #if line_split[0][:5] == "GNG13": #GNG13 comes first in the iteration...
            #for p, c in enumerate(cell_indices_in_datafile):
                    #if float(line_split[c+1]) >= threshold:
                        #cells_in_neuron_cluster[p].insert(1, line_split[0])
        #if line_split[0][:4] == "GNG8":
            #for p, c in enumerate(cell_indices_in_datafile):
                    #if float(line_split[c+1]) >= threshold:
                        #if cells_in_neuron_cluster[p][1] == "GNG13":
                            #cells_in_neuron_cluster[p][1] = "GNG8 + GNG13" #GNG8 comes firse
                        #else:
                            #print(cells_in_neuron_cluster[p][1])
                            #cells_in_neuron_cluster[p].insert(1, line_split[0])

for i in cells_in_neuron_cluster:
    if str(i[1])[:3] != "GNG":
        i.insert(1, "no GNG expression")
    write_to_file(i)
            


# In[2]:


### heatmap
import copy

unique_OR_list = ["VN1R1", "OR10H1", "OR2B11", "OR8D4", "OR8A1", "OR4F6", "OR2F1", "OR2AT4", "OR4N5", "OR4D9", "OR56B1", "OR10Z1", "OR3A3", "OR5AN1", "OR2A42", "OR2A1", "OR6C4", "OR7C1", "OR5AU1", "OR52I1", "OR10G3", "OR11G2", "OR2A25", "OR5A2", "OR2C3", "OR5T1", "OR4F15", "OR13A1", "OR6A2", "OR51M1", "OR7A5", "OR7D4", "OR5V1", "OR8B3", "OR1I1", "OR52N2", "OR6N1", "OR2V2", "OR5B2", "OR8G5", "OR5K1", "OR10A6", "OR2F2", "OR9K2", "OR7A17", "OR2V1", "OR8J3", "OR1L8", "OR1M1", "OR5D14", "OR56A1", "OR2AG2", "OR52I2", "OR6B1", "OR51L1", "OR52K1", "OR52H1", "OR9I1", "OR52N4", "OR1D2", "OR9G1", "OR5A1", "OR3A1", "OR4A16", "OR51E2", "OR6C74", "OR6C6", "OR9Q1", "OR6C2", "OR2M2", "OR2Z1", "OR4K14", "OR2M4", "OR8D1", "OR6C76", "OR52E5", "OR8B2", "OR51G2", "OR10A5", "OR51B5", "OR2L3", "OR2AP1", "OR6C1", "OR7E24", "OR5M10", "OR6N2", "OR1C1", "OR8H2", "OR10G4", "OR2M3", "OR2G6", "OR7D2", "OR5D3P", "OR8G1", "OR5M3", "OR10A2", "OR7G2", "OR6J1", "OR3A2", "OR1E1", "OR56A4", "OR4K13", "OR4M1", "OR2I1P", "OR9Q2", "OR10K1", "OR51E1", "OR7A10", "OR14J1", "OR51V1", "OR9G4", "OR56A3", "OR8B8", "OR6C75", "OR5AS1", "OR8H1", "OR52A1", "OR10G7", "OR5T2", "OR2J3", "OR4D10", "OR10H5", "OR4D1", "OR2T10", "OR11A1", "OR2L2", "OR1A1", "OR52E4", "OR52N5", "OR8K3", "OR4K2", "OR2AJ1", "OR2A5", "OR6F1", "OR11H6", "OR2A12", "OR10K2", "OR4P4", "OR6Y1", "OR8B12", "OR52A5"]
#141
temp = [0] * len(unique_OR_list)
#141
array = []
for i in range(len(unique_OR_list)):
    array.append(copy.copy(temp))
#141x141
    
# make array of OR-co expressions sorted by the unique_OR_list    
for working_OR in unique_OR_list:
    with open("for_heamap.txt") as f:
        for line in f:
            line_split = line.split("\n")[0].split("\t")
            while "" in line_split:
                line_split.pop(line_split.index(""))
            if working_OR in line_split:
                for i in line_split:
                    if i != working_OR:
                        array[unique_OR_list.index(working_OR)][unique_OR_list.index(i)] += 1

list_of_sums = [] # index, sum

for x, i in enumerate(array):
    list_of_sums.append([x, sum(i)])


list_of_sums = sorted(list_of_sums, key=lambda x: x[1], reverse = True)
sorted_array = []
for i in list_of_sums:
    temp = []
    for x, y in enumerate(range(len(list_of_sums))):
        temp.append(array[i[0]][list_of_sums[x][0]])
    sorted_array.append(temp)
    
with open('heatmap.csv', mode='w') as f:
    f_writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for i in sorted_array:
        f_writer.writerow(i)


# In[4]:


### compare gene expression in all clusters
#genes_of_interest = ["VN1R1","GPRC5A" ,"GPRC5B","GPRC5C","GPRC5D","GPRC5D-AS1","ADIPOR1", "GABBR1", "DRD2", "ADGRL3", "TMEM181", "GNG8", "GNG13","GNAL","CNGA2","CNGA4","CNGB1","ADCY3", "GFY", "RTP1", "RTP2", "TUBB3"]
#genes_of_interest = ["TAAR1", "TAAR2", "TAAR5", "TAAR6", "TAAR8", "TAAR9", "VN1R1", "VN1R2", "VN1R3", "VN1R5"]
genes_of_interest = ["VN1R2", "VN1R3", "VN1R5"]


with open(data_file) as f:
    for x, i in enumerate(f):
        line_split = i.split("\n")[0].split("\t")
        if x == 0:
            list_of_cells = line_split
            break

final_output =[]
clusters = []

with open(cell_cluster_correlation_file) as f:
    for row_nr, row in enumerate(f):
        if row_nr > 0:
            line_split = row.split("\n")[0].split("\t")
            if line_split[-1] not in clusters:
                clusters.append(line_split[-1])

for GOI in genes_of_interest:
    output = []
    print(GOI)
    for i in range(len(clusters)):
        cells_counted = [] # list of cells counted to avoid double counting due to multiple receptors
        output.append([clusters[i], 0])  # cluster, count

        cells_in_cluster = []
        with open(cell_cluster_correlation_file) as f:
            for row_nr, row in enumerate(f):
                line_split = row.split("\n")[0].split("\t")
                if line_split[-1] == str(clusters[i]):
                    cells_in_cluster.append(line_split[0])

            output[i].append(len(cells_in_cluster))

            cell_idxs = []
            for cell in cells_in_cluster:
                cell_idxs.append(list_of_cells.index(cell))            

            with open(data_file) as datafile:
                for j in datafile:
                    line_split = j.split("\n")[0].split("\t")
                    #if line_split[0][:2] == "OR":
                    if line_split[0] == GOI:
                        #if line_split[0][2].isdigit() == True:
                        for index in cell_idxs:
                            if float(line_split[index+1]) >= threshold:
                                if index not in cells_counted:
                                    cells_counted.append(index)
                                    output[i][1] +=1

            output[i].append(100*output[i][1]/output[i][2])

    with open(GOI + "_in_all_groups.csv", mode='w') as file:
        file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(["Cluster", "Cells w receptor", "Cells w/o receptor", "Percentage w receptor"])
        for i in output:
            file_writer.writerow(i)
    
    temp = []
    temp.append(GOI)
    for i in output:
        temp.append(output[3])
    final_output.append(temp)

 


# In[5]:


import csv
output_matrix = []
#genes_of_interest = ["VN1R1","GPRC5A" ,"GPRC5B","GPRC5C","GPRC5D","GPRC5D-AS1","ADIPOR1", "GABBR1", "DRD2", "ADGRL3", "TMEM181", "GNG8", "GNG13","GNAL","CNGA2","CNGA4","CNGB1","ADCY3", "GFY", "RTP1", "RTP2", "TUBB3"]
genes_of_interest = ["TAAR1", "TAAR2", "TAAR5", "TAAR6", "TAAR8", "TAAR9", "VN1R1", "VN1R2", "VN1R3", "VN1R5"]


for x, i in enumerate(genes_of_interest):
    with open(i+"_in_all_groups.csv") as datafile:
        for y, j in enumerate(datafile):
            line_split = j.split("\n")[0].split(",")
            if x == 0:
                output_matrix.append([line_split[0]]) #cell cluster label
            if y == 0:
                output_matrix[y].append(i) #gene name
            else:
                output_matrix[y].append(line_split[-1]) # percentage
            
with open("alltogether.csv", mode='w') as file:
    file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for i in output_matrix:
        file_writer.writerow(i)


# In[10]:


classes_OR = []
classes_OR2 = []
with open(data_file) as f:
    for x, i in enumerate(f):
        line_split = i.split("\n")[0].split("\t")
        if line_split[0][:2] == "OR":
            if line_split[0][2].isdigit() == True:
                for p, c in enumerate(cell_indices_in_datafile):
                    if float(line_split[c+1]) >= threshold:
                        temp = line_split[0][2:]
                        for w, q in enumerate(temp):
                            if q.isdigit() == False:
                                classes_OR.append(temp[:w])
                                classes_OR2.append(line_split[0])
                                break
print(classes_OR)
print(classes_OR2)


# In[9]:


print(len(classes_OR))


# In[ ]:




