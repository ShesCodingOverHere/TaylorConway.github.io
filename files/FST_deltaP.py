from __future__ import division
import math

# read in files for p1 and p2
InFile1=open("/nobackup/rogers_research/Taylor/OnceMoreWithFeeling/WithZerosAlleleFreqsIslandDYAK.out",'r')
InFile2=open("/nobackup/rogers_research/Taylor/OnceMoreWithFeeling/WithZerosAlleleFreqsMainlandDYAK.out",'r')
deltaPout=open("/nobackup/rogers_research/Taylor/OnceMoreWithFeeling/DYAKDeltaPWithZeros.out",'w')

p1Dict={}
p2Dict={}

def readDict(Infile, mydict):
    for i in Infile:
        A=i.split()
        span=A[0].split('.')
        chromosome=span[0]
        start=int(span[1])
        stop=int(span[2])
        B=A[1].split('/')
        num=int(B[0])
        denom=int(B[1])
        mykey=A[0]
        if float(B[1])==0:
            p=0
        else:
            p=float(B[0])/float(B[1])
        mydict[mykey]=[chromosome,start,stop,num,denom,p]
readDict(InFile1,p1Dict)
readDict(InFile2,p2Dict)

for key in p1Dict:
    chrom=p1Dict[key][0]
    start=p1Dict[key][1]
    stop=p1Dict[key][2]
    num1=p1Dict[key][3]
    denom1=p1Dict[key][4]
    num2=p2Dict[key][3]
    denom2=p2Dict[key][4]
    freq1=p1Dict[key][5]
    freq2=p2Dict[key][5]

    #midpoint
    mutStart=p1Dict[key][1]
    mutStop=p1Dict[key][2]
    mutMidpoint = str(math.trunc((int(mutStart) + int(mutStop)) / 2))
    #FST
    if (denom1+denom2)!=0:
        p_bar=(num1+num2)/(denom1+denom2)
        q_bar=1-p_bar
        c_1=denom1/(denom1+denom2)
        c_2=denom2/(denom1+denom2)
        p_1=freq1
        p_2=freq2
        q_1=1-p_1
        q_2=1-p_2
        totalnum=num1+num2
        totaldenom=denom1+denom2
        totalfreq=totalnum/totaldenom

        if (p_bar*q_bar)==0:
            fst=0
        else:
            fst=((p_bar*q_bar)-((p_1*q_1*c_1)+(p_2*q_2*c_2)))/(p_bar*q_bar)

        # Delta P
        deltaP=p1Dict[key][5]-p2Dict[key][5]
        deltaPout.write(chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+str(mutMidpoint)+"\t"+str(num1)+"\t"+str(denom1)+"\t"+str(freq1)+"\t"+str(num2)+"\t"+str(denom2)+"\t"+str(freq2)+"\t"+str(deltaP)+"\t"+str(fst)+"\t"+str(totalfreq)+"\n")
