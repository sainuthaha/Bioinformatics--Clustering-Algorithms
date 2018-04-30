import pandas as pd
import numpy as np
import math
import random
from decimal import Decimal
from scipy.spatial import distance

#calculate entropy of a cluster
def entropy(clusterlength,classnumber):
  if clusterlength==0 or classnumber == 0:
    return 0
  tem=0
  tem= float(classnumber)/float(clusterlength)
  prob=tem * math.log(tem,2)
  return prob

#find purity,precision,recall of the cluster
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

def euclideandistance(cent,ele):
  dist = distance.euclidean(cent,ele)
  return dist

data=pd.read_csv('data.csv')
print(data.shape)

k=5
arr=np.array(data)
np.delete(arr,np.s_[0:2],axis=1)
i=0
size= len(arr)


#choosing random initial centroids
num=np.zeros(5)
for i in range(0,5):
  num[i]=random.randint(0,size)
  
col= arr.shape[1]
c=np.zeros((5,col))

for i in range(0,5):
  index=int(num[i])
  c[i]=arr[index]

#calculate the distance of each point from initial centroids and find closest centroid
cluster=np.ones(size)
for i in range(0,size):
	  minimum=euclideandistance(c[0],arr[i])
	  cluster[i]=0
	  for j in range(1,5):
		    dist=euclideandistance(c[j],arr[i])
		    if(dist<minimum):
			       minimum=dist;
			       cluster[i]=j

#Update clusters and find new centroids			    
while(1):
	   flag=0;
	   s0=[]
	   s1=[]
	   s2=[]
	   s3=[]
	   s4=[]
           s0indices=[]
           s1indices=[]
           s2indices=[]
           s3indices=[]
           s4indices=[]

	   for i in range(0,size):
		     clu=cluster[i];
		     if clu== 0:
		       s0.append(arr[i])
                       s0indices.append(i)
		     elif clu==1:
		       s1.append(arr[i])
                       s1indices.append(i)
		     elif clu==2:
		       s2.append(arr[i])
                       s2indices.append(i)
		     elif clu==3:
		       s3.append(arr[i])
                       s3indices.append(i)
                     elif clu==4:
                       s4.append(arr[i])
                       s4indices.append(i)		    

	   s0size=len(s0)
	   s1size=len(s1)
	   s2size=len(s2)
	   s3size=len(s3)
           s4size=len(s4)
    
	   sizearr=np.zeros(5)
	   sizearr[0]=s0size
	   sizearr[1]=s1size
	   sizearr[2]=s2size
	   sizearr[3]=s3size
	   sizearr[4]=s4size
	   sumarr=np.zeros((5,col));
	   
	   k=0  
	   for j in range(0,col):
		s=0
		for i in range(0,s0size):
		   s=s+s0[i][j]
	      
		sumarr[0][j]=s/s0size;
		s=0
		 
		for i in range(0,s1size):
		   s=s+s1[i][j]
		sumarr[1][j]=s/s1size
		s=0
		for i in range(0,s2size):
		   s=s+s2[i][j]
		sumarr[2][j]=s/s2size
		s=0
		for i in range(0,s3size):
		   s=s+s3[i][j]
		sumarr[3][j]=s/s3size
                s=0
	        for i in range(0,s4size):
                   s=s+s4[i][j]
	 
                sumarr[4][j]=s/s4size	  

	   count=0;
	   clus=np.zeros(size)
	   for i in range(0,size):
			    minimum=euclideandistance(sumarr[0],arr[i])
			    clus[i]=0;
			    for j in range(1,5):
                                      
				      dist=euclideandistance(sumarr[j],arr[i])
				                                       
				      if(dist<minimum):
					       count=count+1
					
       
					       minimum=dist;
					       clus[i]=j
	  
          	     
	   for i in range(0,size):
	     if clus[i]!=cluster[i]:
		cluster[i]=clus[i]
		flag=1
	  
	   if flag==0:
	     #print "Final cluster"
	     #print cluster
	     break
	   #else:
	    #print "Current cluster"
	    # print clus
	   
print "Cluster 0:"
print s0indices
print "Cluster 1:"
print s1indices
print "Cluster 2:"
print s2indices
print "Cluster 3:"
print s3indices
print "Cluster 4:"
print s4indices
print "Number of elements in each cluster: ", sizearr
labels=pd.read_csv('labels.csv')
labelarr=np.array(labels)
np.delete(labelarr,np.s_[0:2],axis=1)
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

