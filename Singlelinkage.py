import pandas as pd
import numpy as np
import random
from decimal import Decimal
from scipy.spatial import distance

#Calculation of initial distance matrix
def calc_dist(arr) :
	size= len(arr)
	dist_matrix = np.zeros((size,size))
	for i in range(size-1):
		for j in range(i+1,size):
			dist_matrix[i][j] = dist_matrix[j][i]= euclideandistance(arr[i],arr[j])
	return dist_matrix

#Returns the indices of the clusters with the smallest distance between them
def find_min(dist_matrix) :
	length = dist_matrix.shape[1]
	min_dist = dist_matrix[0][1]
	index_i = 0
	index_j = 1
	for i in range(0,length-1):
		for j in range(i+1,length):
			if(dist_matrix[i][j] < min_dist):
				min_dist = dist_matrix[i][j]
				index_i = i
				index_j = j
	if(index_i < index_j):
		return index_i,index_j
	else:
		return index_j,index_i

#Returns euclidean distance between two data points 
def euclideandistance(cent,ele):
	dist = distance.euclidean(cent,ele)
	return dist

#Returns entropy of a given cluster
def entropy(clusterlength,classnumber):
  if clusterlength==0 or classnumber == 0:
    return 0
  tem=0
  tem= float(classnumber)/float(clusterlength)
  prob=tem * math.log(tem,2)
  return prob

#Returns purity,precision,recall of a cluster
def evalclus(cl,labelarr,l,overallsum):
	clu = np.zeros(5)
	for i in range(0,len(cl)):   
		if labelarr[cl[i]] == 0:
 			clu[0] += 1
		elif labelarr[cl[i]] == 1:
       			clu[1] += 1
  		elif labelarr[cl[i]] == 2:
     			clu[2] += 1
   		elif labelarr[cl[i]] == 3:
       			clu[3] += 1
   		elif labelarr[cl[i]] == 4:
       			clu[4] += 1
	p=0
	for i  in range(5):
		p+= entropy(len(cl),clu[i])
	largest=np.amax(clu)
	overallsum+=largest
	index = np.argmax(clu)
	return float(largest)/float(len(cl)),float(largest)/float(l[index]), -p,overallsum

#Reading the data
data=pd.read_csv('data.csv')
arr=np.array(data)
np.delete(arr,np.s_[0:2],axis=1)
i=0
size= len(arr)
clus = np.arange(size)
dist_matrix = calc_dist(arr)

#Merging clusters until only 5 clusters are remaining
while(dist_matrix.shape[1] > 5):
	i,j = find_min(dist_matrix)
	for x in range(len(clus)):
		if clus[x] == j:
			clus[x] = i
		elif clus[x] > j:
			clus[x] -= 1
	dist_matrix = np.delete(dist_matrix,(j),axis=0)
	dist_matrix = np.delete(dist_matrix,(j),axis=1)
	for p in range(dist_matrix.shape[1]):
		dist_matrix[i][p] = float('Inf')
		dist_matrix[p][i] = float('Inf')
	for x in range(len(clus)):
		if(clus[x] == i) :
			for y in range(len(clus)):
				if(clus[y] != clus[x]):
					dist = euclideandistance(arr[y],arr[x])
					if(dist < dist_matrix[clus[x]][clus[y]]):
						dist_matrix[clus[x]][clus[y]]= dist_matrix[clus[y]][clus[x]] = dist
cl0 = []
cl1 = []
cl2 = []
cl3 = []
cl4 = []
for x in range(size):
	if clus[x] == 0:
		cl0.append(x)
	elif clus[x] == 1:
		cl1.append(x)
	elif clus[x] == 2:
		cl2.append(x)
	elif clus[x] == 3:
		cl3.append(x)
	elif clus[x] == 4:
		cl4.append(x)
print "Cluster0:"
print cl0
print "Cluster1:"
print cl1
print "Cluster2:"
print cl2
print "Cluster3:"
print cl3
print "Cluster4:"
print cl4
print "Number of elements in the clusters: " ,len(cl0), len(cl1), len(cl2), len(cl3), len(cl4)

labels=pd.read_csv('labels.csv')
labelarr=np.array(labels)
labelsize=len(labelarr)

l0=0
l1=0
l2=0
l3=0
l4=0
for i in range(0,labelsize):
  if labelarr[i] == 'BRCA':
   labelarr[i]=0
   l0=l0+1
  elif labelarr[i]== 'LUAD':
   labelarr[i]=1
   l1=l1+1
  elif labelarr[i]== 'PRAD':
   labelarr[i]=2
   l2=l2+1
  elif labelarr[i]== 'KIRC':
   labelarr[i]=3
   l3=l3+1
  elif labelarr[i]== 'COAD':
   labelarr[i]=4
   l4=l4+1

l = [l0,l1,l2,l3,l4]
overallsum=0
pur, recall, entropy0,overallsum = evalclus(cl0,labelarr,l,overallsum)
print "Purity of cluster0: " , pur , " Recall :" , recall 
pur, recall, entropy1,overallsum  = evalclus(cl1,labelarr,l,overallsum)
print "Purity of cluster1: " , pur , " Recall :" , recall 
pur, recall, entropy2,overallsum  = evalclus(cl2,labelarr,l,overallsum)
print "Purity of cluster2: " , pur , " Recall :" , recall 
pur, recall, entropy3,overallsum  = evalclus(cl3,labelarr,l,overallsum)
print "Purity of cluster3: " , pur , " Recall :" , recall 
pur, recall, entropy4,overallsum = evalclus(cl4,labelarr,l,overallsum)
print "Purity of cluster4: " , pur , " Recall :" , recall
print "Overall Purity: ",  float(overallsum)/float(801)
en=float(entropy0 * len(cl0) +entropy1*len(cl1) +entropy2*len(cl2)+entropy3*len(cl3) +entropy4*len(cl4))/float(801)
print "Overall Entropy: ", en
