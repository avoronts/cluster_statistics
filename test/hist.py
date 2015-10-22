#!/usr/bin/python

try: 
   f=open('hist.dat','rb')
   stat = []
   while 1 == 1:
     try:
        i = f.next()
        dat = [j for j in i.replace('(','').replace(')','').replace('+','').split()]
#        print dat
        if dat[2] == '1':
           if dat[3] == dat[4] == '1':
              stat.append(int(dat[1]))
        i = f.next()
     except:
        break
   f.close()

   stat = [j*0.003  for j in stat if j<100]
   f=open('dat.dat','wb')
   f.write(str(stat))
   f.close()
except:
   pass
   
f=open('dat.dat','rb')
time =  [float(j) for j in f.readline().replace('[','').replace(']','').split(',')]
f.close()

try:
   import pylab as pl
   pl.hist(time,30)
   #pl.xlim(0, 200)
   pl.show()
except:
   pass
