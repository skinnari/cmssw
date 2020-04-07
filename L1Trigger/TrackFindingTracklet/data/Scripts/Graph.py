#!/usr/bin/env python

import math

fi = open("memorymodules.dat","r")

memorymodules=[]
memorycounts=[0,0,0,0,0,0,0,0,0,0,0]

for line in fi :
    splitline=line.split(" ")
    #print "Line:",splitline
    lenold=len(memorymodules)
    if (splitline[0]=="InputLink:") :
        memorymodules.append([1,splitline[1]])
    #if (splitline[0]=="StubsByLayer:") :
    #    memorymodules.append([2,splitline[1]])
    #if (splitline[0]=="StubsByDisk:") :
    #    memorymodules.append([2,splitline[1]])
    if (splitline[0]=="AllStubs:" or splitline[0]=="VMStubsTE:" or splitline[0]=="VMStubsME:") :        
        basename=splitline[1].split("n")[0]
        number=0
        if (len(splitline[1].split("n"))>1) :
            number=splitline[1].split("n")[1]
        found=False
        for x in memorymodules :
            if (basename==x[1]) :
                found=True
                if (number>x[2]) :
                    x[2]=number
        if (not found) :
            memorymodules.append([2,basename,number])
    if (splitline[0]=="StubPairs:") :
        memorymodules.append([3,splitline[1]])
    if (splitline[0]=="TrackletParameters:") :
        memorymodules.append([4,splitline[1]])
    if (splitline[0]=="TrackletProjections:") :
        #if "From" not in splitline[1] :
        memorymodules.append([4,splitline[1]])
        #else :
        #    memorymodules.append([6,splitline[1]])
    if (splitline[0]=="AllProj:") :
        memorymodules.append([5,splitline[1]])
    if (splitline[0]=="VMProjections:") :
        memorymodules.append([5,splitline[1]])
    if (splitline[0]=="CandidateMatch:") :
        memorymodules.append([6,splitline[1]])
    if (splitline[0]=="FullMatch:") :
        if "From" not in splitline[1] :
            basename=splitline[1].split("n")[0]
            number=0
            if (len(splitline[1].split("n"))>1) :
                number=splitline[1].split("n")[1]
            found=False
            for x in memorymodules :
                if (basename==x[1]) :
                    found=True
                    if (number>x[2]) :
                        x[2]=number
            if (not found) :
                memorymodules.append([7,basename,number])
    if (splitline[0]=="TrackFit:") :
        memorymodules.append([8,splitline[1]])
    if (splitline[0]=="CleanTrack:") :
        memorymodules.append([9,splitline[1]])

    if (lenold != len(memorymodules)) :    
        memorycounts[memorymodules[len(memorymodules)-1][0]-1]+=1

#for x in memorymodules :
#    print x

print "Memorycounts:",memorycounts



fi = open("processingmodules.dat","r")

processingmodules=[]
processingcounts=[0,0,0,0,0,0,0,0,0,0]

for line in fi :
    splitline=line.split(" ")
    #print "Line:",splitline
    lenold=len(processingmodules)
    if (splitline[0]=="VMRouter:") :
        processingmodules.append([1,splitline[1]])
    if (splitline[0]=="TrackletEngine:") :
        processingmodules.append([2,splitline[1]])
    if (splitline[0]=="TrackletCalculator:") :
        processingmodules.append([3,splitline[1]])
    if (splitline[0]=="ProjectionRouter:") :
        processingmodules.append([4,splitline[1]])
    if (splitline[0]=="MatchEngine:") :
        processingmodules.append([5,splitline[1]])
    if (splitline[0]=="MatchCalculator:") :
        processingmodules.append([6,splitline[1]])
    if (splitline[0]=="DiskMatchCalculator:") :
        processingmodules.append([6,splitline[1]])
    if (splitline[0]=="FitTrack:") :
        processingmodules.append([7,splitline[1]])
    if (splitline[0]=="PurgeDuplicate:") :
        processingmodules.append([8,splitline[1]])

    if (lenold != len(processingmodules)) :    
        processingcounts[processingmodules[len(processingmodules)-1][0]-1]+=1

#for x in processingmodules :
#    print "processingmodules : ",x

print "Processing counts:",processingcounts

#           IL    AS   SP  TPAR  AP   CM    FM   TF   CT
xmemories=[0.005,0.16, 0.31,0.47,0.60,0.72, 0.80, 0.91,0.97]
dxmemories=[0.07,0.045,0.055,0.055,0.05,0.035,0.045,0.025,0.025]
memories=[0,0,0,0,0,0,0,0,0]

fo = open("diagram.dat","w")

for module in memorymodules :
    num=memories[module[0]-1]
    memories[module[0]-1]+=1
    y=(0.5+num)/(memorycounts[module[0]-1])
    x=xmemories[module[0]-1]
    dx=dxmemories[module[0]-1]
    module.append([x,y,x+dx,y])
    fo.write("Memory "+module[1])
    fo.write(" "+str(x)+" "+str(y)+" "+str(x+dx)+" "+str(y)+"\n")
    #print module


#            VMR    TE   TC    PR    ME    MC    FT   PD
xprocessing=[0.095, 0.24,0.40, 0.55, 0.67, 0.76, 0.87,0.94,0.98]
dxprocessing=[0.035,0.06,0.035,0.035,0.035,0.035,0.03,0.02]
processing=[0,0,0,0,0,0,0,0]

for module in processingmodules :
    num=processing[module[0]-1]
    processing[module[0]-1]+=1
    y=(0.5+num)/(processingcounts[module[0]-1])
    x=xprocessing[module[0]-1]
    dx=dxprocessing[module[0]-1]
    module.append([x,y,x+dx,y])
    fo.write("Process "+module[1].split("\n")[0])
    fo.write(" "+str(x)+" "+str(y)+" "+str(x+dx)+" "+str(y)+"\n")
    #print module

fi = open("wires.dat","r")
for line in fi :
    memory=line.split(" ")[0]
    if (memory[len(memory)-2]=="n" or memory[len(memory)-3]=="n") :
        memory=memory.split("n")[0]
    inprocess=line.split("input=> ")[1].split(" ")[0].split(".")[0]
    #print "line:",line
    outprocess=""
    if len(line.split("output=> "))>1 :
        outprocess=line.split("output=> ")[1].split("\n")[0].split(".")[0]

    #print memory+" in: "+inprocess+" out: "+outprocess

    #Find the memory module
    memmodule=[]
    for mem in memorymodules :
        if (mem[1]==memory) :
            memmodule=mem
    if (memmodule==[]) :
        print "Could not find memorymodule: ",memory
        print "Will terminate"
        exit()

    last=len(memmodule)-1
        
    #Find the input processing module    
    if (inprocess!="") :
        procmodule=[]        
        for proc in processingmodules :
            #print "Comparing ",proc[1].split("\n")[0]," - ",inprocess
            if (proc[1].split("\n")[0]==inprocess) :
                procmodule=proc
        if (procmodule==[]) :
            print "Could not find in processingmodule: ",inprocess
            print "Will terminate"
            exit()
        if (memmodule[0]==procmodule[0]+1 or memmodule[0]==procmodule[0]) :
            fo.write("Line "+str(procmodule[2][2])+" "+str(procmodule[2][3]))
            fo.write(" "+str(memmodule[last][0])+" "+str(memmodule[last][1])+"\n")
        

    #Find the output processing module    
    if (outprocess!="") :    
        procmodule=[]        
        for proc in processingmodules :
            #print "Comparing ",proc[1].split("\n")[0]," - ",outprocess
            if (proc[1].split("\n")[0]==outprocess) :
                procmodule=proc
        if (procmodule==[]) :
            print "Could not find outprocessingmodule: ",outprocess
            print "Will terminate"
            exit()
        #if (memmodule[0]==procmodule[0] or memmodule[0]==procmodule[0]-1) :
        if (memmodule[0]==procmodule[0]) :
            fo.write("Line "+str(procmodule[2][0])+" "+str(procmodule[2][1]))
            fo.write(" "+str(memmodule[last][2])+" "+str(memmodule[last][3])+"\n")
