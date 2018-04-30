import pandas as pd
import numpy as np
import random
from decimal import Decimal
from scipy.spatial import distance
import math

#Finds euclidean distance between 2 data points
def euclideandistance(cent,ele):
   dist = distance.euclidean(cent,ele)
   return dist

#Creates labels for each of the n data points
def classlbl(start, end):
	labels = []
	for i in range(start, end+1):
		li=[]
                li.append(i)
		labels.append(li)
	return labels

#Calculates entropy for a cluster
def entropy(clusterlength,classnumber):
  if clusterlength==0 or classnumber == 0:
    return 0
  tem=0
  tem= float(classnumber)/float(clusterlength)
  prob=tem * math.log(tem,2)
  return prob

#Calculates purity,precision,recall for a cluster
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

#Read the data set
data=pd.read_csv('data.csv')
arr=np.array(data)
size= len(arr)

m=[]           #
m_labels=[]

for i in range(0,size):
	  li=[]
	  for j in range(0,i):
		   d=euclideandistance(arr[j],arr[i])
		   li.append(d)
	  m.append(li)
m_labels=classlbl(0,800)

#Returns the indices of the clusters having minimum distance between them
def min_val(dist_mat):
	min_cell = float("inf")
	x, y = -1, -1
	for i in range(len(dist_mat)):
		for j in range(len(dist_mat[i])):
			if dist_mat[i][j] < min_cell:
				min_cell = dist_mat[i][j]
				x, y = i, j
	return x, y

#merges the labels of clusters having minimum distance
def mergeclus(labels, a, b):
	if b < a:
		a, b = b, a
        labels[a] += labels[b]
        del labels[b]

#Updates the distance matrix by finding the distance of the merged cluster from the other clusters
def update_dist(dist_mat, a, b):
	if b < a:
		a, b = b, a
	row = []
	for i in range(0, a):
                avg=float(dist_mat[a][i]+dist_mat[b][i])/float(2)
		row.append(avg)
	dist_mat[a] = row
	for i in range(a+1, b):
		dist_mat[i][a] = float(dist_mat[i][a]+dist_mat[b][i])/float(2)
	
	for i in range(b+1, len(dist_mat)):
		dist_mat[i][a] = float(dist_mat[i][a]+dist_mat[i][b])/float(2)
		del dist_mat[i][b]
	
	del dist_mat[b]

#Returns a list containing clusters obtained by using upgma
def UPGMA(dist_mat, labels):
	while len(labels) > 5:
		x, y = min_val(dist_mat)
		update_dist(dist_mat, x, y)
		mergeclus(labels, x, y)
	return labels

labs=[]
labs=UPGMA(m,m_labels)
print "Cluster0:"
print labs[0]
print "Cluster1:"
print labs[1]
print "Cluster2:"
print labs[2]
print "Cluster3:"
print labs[3]
print "Cluster4:"
print labs[4]

l0=len(labs[0])
l1=len(labs[1])
l2=len(labs[2])
l3=len(labs[3])
l4=len(labs[4])

s0indices=np.asarray(labs[0])
s1indices=np.asarray(labs[1])
s2indices=np.asarray(labs[2])
s3indices=np.asarray(labs[3])
s4indices=np.asarray(labs[4]) 
print "Size of clusters: ", (l0,l1,l2,l3,l4)

labels=pd.read_csv('labels.csv')
labelarr=np.array(labels)
np.delete(labelarr,np.s_[0:2],axis=1)
labelsize=len(labelarr)

a0=0
a1=0
a2=0
a3=0
a4=0
for i in range(0,labelsize):
  if labelarr[i] == 'BRCA':
   labelarr[i]=0
   a0=a0+1
  elif labelarr[i]== 'LUAD':
   labelarr[i]=1
   a1=a1+1
  elif labelarr[i]== 'PRAD':
   labelarr[i]=2
   a2=a2+1
  elif labelarr[i]== 'KIRC':
   labelarr[i]=3
   a3=a3+1
  elif labelarr[i]== 'COAD':
   labelarr[i]=4
   a4=a4+1

l = [a0,a1,a2,a3,a4]
overallsum=0
pur, recall, entropy0,overallsum = evalclus(s0indices,labelarr,l,overallsum)
print "Cluster0- Purity: " , pur , "Precision:" , pur, "Recall :" , recall , "Entropy: " , entropy0
pur, recall, entropy1,overallsum  = evalclus(s1indices,labelarr,l,overallsum)
print "Cluster1- Purity: " , pur , "Precision:" , pur, "Recall :" , recall , "Entropy: " , entropy1
pur, recall, entropy2,overallsum  = evalclus(s2indices,labelarr,l,overallsum)
print "Cluster2- Purity: " , pur , "Precision:" , pur, "Recall :" , recall , "Entropy: " , entropy2
pur, recall, entropy3,overallsum  = evalclus(s3indices,labelarr,l,overallsum)
print "Cluster3- Purity: " , pur , "Precision:" , pur, "Recall :" , recall , "Entropy: " , entropy3
pur, recall, entropy4,overallsum = evalclus(s4indices,labelarr,l,overallsum)
print "Cluster4- Purity: " , pur , "Precision:" , pur, "Recall :" , recall , "Entropy: " , entropy4
print "Overall Purity: ",  float(overallsum)/float(801)
en=float(entropy0 * len(s0indices) +entropy1*len(s1indices) +entropy2*len(s2indices)+entropy3*len(s3indices) +entropy4*len(s4indices))/float(801)
print "Overall Entropy: ", en

