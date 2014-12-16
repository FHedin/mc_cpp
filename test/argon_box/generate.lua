#!/usr/bin/env lua

-- This script generates the coordinate and forcefield files for a cubic box of argon atoms
-- By modifying the two following variables you can choose the number of particles 
-- and the box size. Choose values carefully because combined with the temperature it defines
-- the phase (solid,liquid,gas) simulated.

nparticl=512
boxa=30

gridl=math.pow(nparticl,1.0/3.0)
shift=boxa/gridl

resid=1
resname="ARR"
atname="AR"
segname="ALL"
segid=1
weight=0.00000

atid=1
x={}
y={}
z={}

for i=1,nparticl
do
    x[i]=0.0
    y[i]=0.0
    z[i]=0.0
end

l=1
for i=1,boxa,shift
do
    for j=1,boxa,shift
    do
        for k=1,boxa,shift
        do
            x[l]=i
            y[l]=j
            z[l]=k
            l=l+1
        end
    end
end

cm={}
cm[1]=0.0
cm[2]=0.0
cm[3]=0.0
for i=1,nparticl
do
    cm[1]=cm[1]+x[i]
    cm[2]=cm[2]+y[i]
    cm[3]=cm[3]+z[i]
end
cm[1]=cm[1]/nparticl
cm[2]=cm[2]/nparticl
cm[3]=cm[3]/nparticl
for i=1,nparticl
do
    x[i]=x[i]-cm[1]
    y[i]=y[i]-cm[2]
    z[i]=z[i]-cm[3]
end

outfile=io.open("ar.cor","w")
str=string.format("*\n*\n*\n%5d\n",nparticl)
outfile:write(str)
for i=1,nparticl
do
    atid=i
    str=string.format("%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",atid,resid,resname,atname,x[i],y[i],z[i],segname,segid,weight)
    outfile:write(str)
end
outfile:close()

attype=1
charge=0.0
mass=39.948000
vdwtype=0
epsi=0.238000
sig=1.705000
beta=0.0
eps14=0.0
sig14=0.0
beta14=0.0

outfile=io.open("forfield.dat","w")

str=string.format("Atoms\t%d\n",nparticl)
outfile:write(str)

for i=1,nparticl
do
    atid=i
    str=string.format("%d\t%s\t%d\t% #f\t% #f\t%d\n",atid,atname,attype,charge,mass,0.0)
    outfile:write(str)
end

outfile:write(string.format("Bonds\t%d\n",0))
outfile:write(string.format("Angles\t%d\n",0))
outfile:write(string.format("Urey-Bradley\t%d\n",0))
outfile:write(string.format("Dihedrals\t%d\n",0))
outfile:write(string.format("Impropers\t%d\n",0))
outfile:write(string.format("Constraints\t%d\n",0))

str=string.format("vdw\t%d\n",nparticl)
outfile:write(str)

for i=1,nparticl
do
    atid=i
    str=string.format("%d\t%d\t% #f\t% #f\t% #f\t% #f\t% #f\t% #f\n",atid,vdwtype,epsi,sig,beta,eps14,sig14,beta14)
    outfile:write(str)
end
outfile:write(string.format("end\n"))

outfile:close()


