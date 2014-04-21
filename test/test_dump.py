#!/usr/bin/python

import os
try:
   home_dir=os.environ['HOME']
except:
   home_dir='e:'

import sys
sys.path += [home_dir+'/bin/model/']
sys.path += [home_dir+'/bin/utils/']
sys.path += [os.getcwd()]
############### end of common block for all programms #######################
#def test()

f=open('dump_clust.40000050','rb')
ff=f.readlines()
f.close()
dd=[]
for i in ff[9:]:
  try:
    dd.append([int(j) for j in i.split()])
  except:
    break

dat={}
clust=[0]+[0 for i in dd]
size=[0]+[0 for i in dd]
gen=[0]+[0 for i in dd]

for i in dd:
#  print i
#  clust[i[0]]=i[2]
#  size[i[0]]=i[3]
#  gen[i[0]]=i[4]
  dat[i[2]] = dat.get(i[2],[])+[i[0]]

print len(dat), len(dd) #, sum([abs(i) for i in size])
stop
for i in sorted(dat.keys()):
   for j in dat[i]:
     if len(dat[i]) <> size[j]:
       print  i, len(dat[i]), size[j] #[size[k] for k in dat[i]]

#  if  len(dat[i])>2:
#    print  i, len(dat[i]), [(size[j],gen[j]) for j in dat[i] if size[j]<>0], sum([size[j] for j in dat[i] if size[j]<>0])

print len(dat), len(dd), sum([abs(i) for i in size])

#for i in dat[10548]:
#  print i,size[i],clust[i], gen[i]

stop
