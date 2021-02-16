# Input polymorhism data for narnavirus, re-direct output.
import pandas as pd
import numpy as np
import math
import sys
file=open('Robin.txt','r+')
data = file.readlines()
sys.stdout = open('Robin.op','wt')
print('Robin data')
print(' ')
#
# Define some procedures.
#
# Define procedure to flip sequence to its complement.
def complement (seq):
    n=len(seq)
    seqcopy=[]
    for i in range (n):
        seqcopy.append('x')
    for i in range (n):
        if seq[i]=='A': seqcopy[i]='T'
        if seq[i]=='T': seqcopy[i]='A'
        if seq[i]=='C': seqcopy[i]='G'
        if seq[i]=='G': seqcopy[i]='C'
    for i in range (n):
        seq[n-i-1]=seqcopy[i]
    return seq
# Function to assign numerical values to codons: base 4, (A,C,G,T)->(0,1,2,3).
def codonval (codon):
    nv=np.zeros(3,dtype=int)
    for i in range (3):
        symb=codon[i]
        if symb=='A': nv[i]=0
        if symb=='C': nv[i]=1
        if symb=='G': nv[i]=2
        if symb=='T': nv[i]=3
    codonval=nv[0]+4*nv[1]+16*nv[2]
    return codonval
# Function to return index of reverse-complement of a codon.
def flip (ind):
    nv=np.zeros(3,dtype=int)
    nvf=np.zeros(3,dtype=int)
    nv[2]=int(ind/16)
    nv[1]=int((ind-16*nv[2])/4)
    nv[0]=ind-16*nv[2]-4*nv[1]
    for i in range (3):
        if nv[i]==0: nvf[i]=3
        if nv[i]==1: nvf[i]=2
        if nv[i]==2: nvf[i]=1
        if nv[i]==3: nvf[i]=0
    flip=nvf[2]+4*nvf[1]+16*nvf[0]
    return flip
# Function to return aminoacid code letter from codon index: use 'Z' for stop codons.   
def gencode (ind):
    if ind in [6,22,38,54]: aa='A'
    if ind in [27,59]: aa='C'
    if ind in [18,50]: aa='D'
    if ind in [2,34]: aa='E'
    if ind in [31,63]: aa='F'
    if ind in [10,26,42,58]: aa='G'
    if ind in [17,49]: aa='H'
    if ind in [27,59]: aa='C'
    if ind in [12,28,60]: aa='I'
    if ind in [0,32]: aa='K'
    if ind in [13,15,29,45,47,61]: aa='L'
    if ind in [44]: aa='M'
    if ind in [16,48]: aa='N'
    if ind in [5,21,37,53]: aa='P'
    if ind in [1,33]: aa='Q'
    if ind in [8,9,25,40,41,57]: aa='R'
    if ind in [7,23,24,39,55,56]: aa='S'
    if ind in [4,20,36,52]: aa='T'
    if ind in [14,30,46,62]: aa='V'
    if ind in [43]: aa='W'
    if ind in [19,51]: aa='Y'
    if ind in [3,11,35]: aa='Z'
    gencode=aa
    return gencode
# Function to count nucleotide changes in mutation from i1 to i2.
def ntchanges (i1,i2):
    nv1=np.zeros(3,dtype=int)
    nv2=np.zeros(3,dtype=int)
    nv1[2]=int(i1/16)
    nv1[1]=int((i1-16*nv1[2])/4)
    nv1[0]=i1-16*nv1[2]-4*nv1[1]
    nv2[2]=int(i2/16)
    nv2[1]=int((i2-16*nv2[2])/4)
    nv2[0]=i2-16*nv2[2]-4*nv2[1]
    nd=0
    for i in range (3):
        if nv1[i]!=nv2[i]: nd=nd+1
    ntchanges=nd
    return ntchanges
# Function to return numbers of ns and s transitions (inn, isn) and ns and s tranversions (inv, isv) of a codon.
def codon_mutns (ind):
    inn=0
    isn=0
    inv=0
    isv=0
    aa=gencode(ind)
    n=np.zeros(3,dtype=int)
    nv=np.zeros(3,dtype=int)
    n[2]=int(ind/16)
    n[1]=int((ind-16*n[2])/4)
    n[0]=ind-16*n[2]-4*n[1]
    for j in range (3):
        for k in range (4):
            nv[0]=n[0]
            nv[1]=n[1]
            nv[2]=n[2]
            nv[j]=k
            ind1=nv[0]+4*nv[1]+16*nv[2]
            aa1=gencode(ind1)
            if nv[j] != n[j]:
                if abs(nv[j]-n[j]) == 2:
                    if aa1!=aa: inn=inn+1
                    if aa1==aa: isn=isn+1
                if abs(nv[j]-n[j]) != 2:
                    if aa1!=aa: inv=inv+1
                    if aa1==aa: isv=isv+1
    return inn,inv,isn,isv
# Function to return numbers of double-syn transitions/transversions (idn, idv) of a codon (excluding stops).
def dble_syns (ind):
    idn=0
    idv=0
    aa=gencode(ind)
    aac=gencode(flip(ind))
    n=np.zeros(3,dtype=int)
    nv=np.zeros(3,dtype=int)
    n[2]=int(ind/16)
    n[1]=int((ind-16*n[2])/4)
    n[0]=ind-16*n[2]-4*n[1]
    for j in range (3):
        for k in range (4):
            nv[0]=n[0]
            nv[1]=n[1]
            nv[2]=n[2]
            nv[j]=k
            ind1=nv[0]+4*nv[1]+16*nv[2]
            aa1=gencode(ind1)
            aa1c=gencode(flip(ind1))
            if nv[j] != n[j]:
                if abs(nv[j]-n[j]) == 2:
                    if aa1==aa and aa1c==aac and aac!='Z': idn=idn+1
                if abs(nv[j]-n[j]) != 2:
                    if aa1==aa and aa1c==aac and aac!='Z': idv=idv+1
    return idn,idv
#
# Main code resumes here.
#
# Determine number of polymorphs npol, lengths of sequences ns[i], 
# maximum length sequence has index im, nm entries.
npol=len(data)
ns=np.zeros((npol),dtype=int)
for i in range (npol):
    ns[i]=len(data[i])-1
nm=max(ns)
im=np.where(ns==nm)[0][0]
print('M=',npol,'3N=',nm)
# Convert input data from strings to lists: seqset[i] is list of bases for polymorph i.
seqset=[]
for i in range (npol):
   seq=[]
   seq.extend(data[i])
   seq.remove('\n')
   seqset.append(seq)
#
# Create list of variant sets of bases at each nucleotide locus.
#
varset_bases=[set() for i in range(nm)]
for i in range (nm):
    for j in range (npol):
        base=seqset[j][i]
        varset_bases[i].add(base)
# For each locus, count occurences of bases. 
base_counts=[]
for i in range (nm):
    base_count=np.zeros(4,dtype=int)
    for j in range (npol):
        base=seqset[j][i]
        if base=='A': base_count[0]=base_count[0]+1
        if base=='C': base_count[1]=base_count[1]+1
        if base=='G': base_count[2]=base_count[2]+1
        if base=='T': base_count[3]=base_count[3]+1
    base_counts.append(base_count)
# Determine reference bases, mutation counts, transition counts.
refbases=[]
mutn_counts=[]
trns_counts=[]
ia=np.zeros(4,dtype=int)
for i in range (nm):
    ia=base_counts[i]
    mx=max(ia)
    ix=list(ia).index(mx)
    if ix==0: refbases.append('A')
    if ix==1: refbases.append('C')
    if ix==2: refbases.append('G')
    if ix==3: refbases.append('T')
    mutn_counts.append(0)
    trns_counts.append(0)
    for j in range (4):
        if j!=ix:
            if base_counts[i][j] != 0: mutn_counts[i]=mutn_counts[i]+1
    if ix==0 and base_counts[i][2]!=0: trns_counts[i]=1
    if ix==1 and base_counts[i][3]!=0: trns_counts[i]=1
    if ix==2 and base_counts[i][0]!=0: trns_counts[i]=1
    if ix==3 and base_counts[i][1]!=0: trns_counts[i]=1
# Print aligned nucleotide sequences, if desired:
#for i in range (nm):
#    print(i,refbases[i],'    ',
#    seqset[0][i],seqset[1][i],seqset[2][i],seqset[3][i],seqset[4][i],
#    seqset[5][i],seqset[6][i],seqset[7][i],seqset[8][i],seqset[9][i],
#    seqset[10][i],seqset[11][i],seqset[12][i],seqset[13][i],seqset[14][i],
#    seqset[15][i],seqset[16][i],seqset[17][i],seqset[18][i],seqset[19][i],
#    seqset[20][i],seqset[21][i],seqset[22][i],seqset[23][i],seqset[24][i],
#    seqset[25][i],seqset[26][i],seqset[27][i],seqset[28][i],seqset[29][i],
#    seqset[30][i],seqset[31][i],seqset[32][i],seqset[33][i],seqset[34][i],
#    seqset[35][i],seqset[36][i],seqset[37][i],seqset[38][i],seqset[39][i],
#    seqset[40][i],seqset[41][i],seqset[42][i],seqset[43][i],seqset[44][i],
#    seqset[45][i],varset_bases[i],base_counts[i])
# Estimate alpha: remember there are twice as many routes to transversions.
# Print counts at each locus if desired.
#for i in range (nm-1):
#    print(i,base_counts[i][0],base_counts[i][1],base_counts[i][2],base_counts[i][3],mutn_counts[i],trns_counts[i],fbvals[i])
ntrns=0
ntrvs=0
for i in range (nm):
    ntrns=ntrns+trns_counts[i]
    ntrvs=ntrvs+mutn_counts[i]
ntrvs=ntrvs-ntrns
alpha=2*ntrns/ntrvs
print(' ')
print('nucleotide-level mutation data')
print('nbases=',nm,'ntrns=',ntrns,'ntrvs=',ntrvs,'alpha=',alpha)
# Look at mutation frequency as a function of position in codon (in forward direction).
n1=0
n2=0
n3=0
nc=int(nm/3)
for i in range (nc):
    n1=n1+mutn_counts[3*i]
    n2=n2+mutn_counts[3*i+1]
    n3=n3+mutn_counts[3*i+2]
print('mutation counts by codon position:',n1,n2,n3)
mutn_rate=(n1+n2+n3)/(nm*npol)
print('mutation rate=',mutn_rate)
# 
# Now switch to reverse complement sequence, if required (comment out for forward sequence):
#
#for i in range (npol):
#    complement (seqset[i])
#complement(refbases)
#print('complementary sequence')
print('forward sequence')
# Now identify codons: assign numbers in range 0-63 using base 4, with weights A=0, C=1, G=2, T=3.
# frameshift=0 for forward, frameshift=2 for complement.
frameshift=0
print('frameshift=',frameshift,'(fwd.=0, comp.=0)')
codonvaluesset=[]
for i in range (npol):
    codonvalues=[]
    for j in range (nm):
        codon=['N','N','N']
        for k in range (3):
            index=frameshift+3*j+k
            if index in range (nm):
                symb=seqset[i][index]
                if symb in ('A','C','G','T'):
                    codon[k]=symb
        if 'N' not in codon:
            codonvalues.append(codonval(codon))
        if 'N' in codon:
            codonvalues.append('X')
    codonvaluesset.append(codonvalues)
# Reference codons.
refcodons=[]
for i in range (nm):
    codon=['N','N','N']
    for j in range (3):
        index=frameshift+3*i+j
        if index in range (nm):
            symb=refbases[index]
            if symb in ('A','C','G','T'):
                codon[j]=symb
    if 'N' not in codon:
        refcodons.append(codonval(codon))
    if 'N' in codon:
        refcodons.append('X')
codonvaluesset.append(codonvalues)
# Print codon assignments, if desired.
#for i in range (nc):
#    print(i,refcodons[i],'   ',
#    codonvaluesset[0][i],codonvaluesset[1][i],codonvaluesset[2][i],codonvaluesset[3][i],
#    codonvaluesset[4][i],codonvaluesset[5][i],codonvaluesset[6][i],codonvaluesset[7][i],
#    codonvaluesset[8][i],codonvaluesset[9][i],codonvaluesset[10][i],codonvaluesset[11][i],
#    codonvaluesset[12][i],codonvaluesset[13][i],codonvaluesset[14][i],codonvaluesset[15][i],
#    codonvaluesset[16][i],codonvaluesset[17][i],codonvaluesset[18][i],codonvaluesset[19][i],
#    codonvaluesset[20][i],codonvaluesset[21][i],codonvaluesset[22][i],codonvaluesset[23][i],
#    codonvaluesset[24][i],codonvaluesset[25][i],codonvaluesset[26][i],codonvaluesset[27][i],
#    codonvaluesset[28][i],codonvaluesset[29][i],codonvaluesset[30][i],codonvaluesset[31][i],
#    codonvaluesset[32][i],codonvaluesset[33][i],codonvaluesset[34][i],codonvaluesset[35][i],
#    codonvaluesset[36][i],codonvaluesset[37][i],codonvaluesset[38][i],codonvaluesset[39][i],
#    codonvaluesset[40][i],codonvaluesset[41][i],codonvaluesset[42][i],codonvaluesset[43][i],
#    codonvaluesset[44][i],codonvaluesset[45][i])
#
# Find stops:
#    forward stops are [TAA,TAG,TGA]=[3,35,11]
#    reverse stops are [TTA,CTA,TCA]=[15,13,7]
#
stopset=[3,35,11]
revstopset=[15,13,7]
print('forward stops')
for i in range (nc):
   nstp=0
   for j in range (npol):
       if codonvaluesset[j][i] in stopset:
           nstp=nstp+1
   if nstp>0:
       print('codon=',i,'nstp=',nstp)
print('comp. stops')
for i in range (nc):
   nstp=0
   for j in range (npol):
       if codonvaluesset[j][i] in revstopset:
           nstp=nstp+1
   if nstp>0:
       print('codon=',i,'nstp=',nstp)
#
# Create list of variant sets of codons at each locus.
#
varset_codons=[set() for i in range(nc)]
for i in range (nc):
    for j in range (npol):
        index=codonvaluesset[j][i]
        varset_codons[i].add(index)
# For each locus, count occurences of codons. 
codon_counts=[]
for i in range (nc):
    codon_count=np.zeros(64,dtype=int)
    for j in range (npol):
        ix=codonvaluesset[j][i]
        if ix in range (64):
            codon_count[ix]=codon_count[ix]+1
    codon_counts.append(codon_count)
# Determine fraction of differences from reference codon.
fvals=[]
for i in range (nc):
    nt=0
    for j in range (npol):
        if codonvaluesset[j][i] in range (64): nt=nt+1
    ref=refcodons[i]
    mx=codon_counts[i][ref]
    fv=(nt-mx)/nt
    fvals.append(fv)
# Count numbers of non-syn, syn, and non-SNP mutations.
N_ns=0
N_s=0
N_ns2=0
N_s2=0
N_ns3=0
N_s3=0
N_tot=0
for i in range (nc):
    ind1=refcodons[i]
    if ind1 in range (64):
        aa1=gencode(ind1)
        for ind2 in varset_codons[i]:
            if ind2 in range (64):    
                if ind2 !=ind1:
                    N_tot=N_tot+1
                    aa2=gencode(ind2)
                    n_diff=ntchanges(ind1,ind2)
                    if n_diff==1:
                        if aa2==aa1: N_s=N_s+1
                        if aa2!=aa1: N_ns=N_ns+1
                    if n_diff==2:
                        if aa2==aa1: N_s2=N_s2+1
                        if aa2!=aa1: N_ns2=N_ns2+1
                    if n_diff==3:
                        if aa2==aa1: N_s3=N_s3+1
                        if aa2!=aa1: N_ns3=N_ns3+1
N_mult=N_s2+N_s3+N_ns2+N_ns3
r=N_ns/N_s
print(' ')
print('codon-level mutations')
print('N_tot=',N_tot,'N_ns=',N_ns,'N_s=',N_s,'N_mult=',N_mult,'N/S ratio=',r)
# Determine expected N_ns/N_s ratio.
num=0
den=0
n=0
for i in range (nc):
    ind=refcodons[i]
    inn,inv,isn,isv=codon_mutns(ind)
    n=n+1
    num=num+alpha*inn+inv
    den=den+alpha*isn+isv
r_th=num/den
N_s_exp=N_tot/(1+r_th)
print('r_th=',r_th,'N_s_exp=',N_s_exp)
print(' ')
# Determine variant counts and single and double synonym counts for double-syn loci.
# List of double-syn codons.
dblsynlist=[1,8,9,33,37,39,40,41,45,47,53,55]
# Compute statistics for activity test:
varcounts=[]
for i in range (nc):
    varcounts.append(len(varset_codons[i])) 
ndbl=0
nndbl=0
ndav=0
nndav=0
fdav=0
fndav=0
for i in range (nc):
    ref=refcodons[i]
    if ref in dblsynlist:
        ndbl=ndbl+1
        ndav=ndav+varcounts[i]-1
        fdav=fdav+fvals[i]
    if ref not in dblsynlist:
        nndbl=nndbl+1
        nndav=nndav+varcounts[i]-1
        fndav=fndav+fvals[i]
ndav=ndav/ndbl
nndav=nndav/nndbl
fdav=fdav/ndbl
fndav=fndav/nndbl
nratio=ndav/nndav
fratio=fdav/fndav
print(' ')
print('Hotspots test')
print('N_d=',ndbl,'N-N_d=',nndbl)
print('<n_d>=',ndav,'<n_nd>=',nndav,'<f_d>=',fdav,'<f_nd>=',fndav)
print('nratio=',nratio,'fratio=',fratio)
print(' ')
# Compute statistics for codon variability test:
print('Mutation frequency test')
print('alpha=',alpha)
nstot=0
ndtot=0
ndbl=0
num=0
den=0
for i in range (nc):
    ref=refcodons[i]
    if ref in range (64):
        aa=gencode(ref)
        aac=gencode(flip(ref))
        inn,inv,isn,isv=codon_mutns(ref)
        idn,idv=dble_syns(ref)
    if ref in dblsynlist and varcounts[i] > 1:
        ndbl=ndbl+1
        for j in varset_codons[i]:
            if (j in range (64) and j!=ref):
                aa1=gencode(j)
                aa1c=gencode(flip(j))
                if aa1==aa: nstot=nstot+1
                if (aa1==aa and aa1c==aac): ndtot=ndtot+1
        num=num+alpha*isn+isv
        den=den+alpha*idn+idv
r_th=num/den
print('N_a=',ndbl,'N_s=',nstot,'N_d=',ndtot,'R=',nstot/ndtot,'R_th=',r_th)
sigma=(r_th*ndtot-nstot)/math.sqrt(nstot)
print('sigma=',sigma)
print(' ')






