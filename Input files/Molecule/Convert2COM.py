import numpy

FileName='TIP5PE.xyz'

with open(FileName,'r') as f:
    data=f.readlines()
    N=int(data[3].strip())
    Element=numpy.empty(N,dtype=str)
    Coordinate=numpy.empty((3,N))
    mass=numpy.empty(N)
    COM=numpy.zeros(3)
    for i in range(N):
        temp=data[5+i].split()
        Element[i]=temp[0].strip()
        if(Element[i]=='C'):
            mass[i]=12.011
        elif(Element[i]=='H'):
            mass[i]=1.008
        elif(Element[i]=='O'):
            mass[i]=15.999
        else:
            mass[i]=0.0
        for j in range(3):
            Coordinate[j,i]=float(temp[j+1].strip())
            COM[j]=COM[j]+mass[i]*Coordinate[j,i]

COM=COM/sum(mass)
for i in range(N):
    Coordinate[:,i]=Coordinate[:,i]-COM

with open('1.txt','w') as f:
    for i in range(N):
        print("%2s%20.15f%20.15f%20.15f"%(Element[i],Coordinate[0,i],Coordinate[1,i],Coordinate[2,i]),file=f)