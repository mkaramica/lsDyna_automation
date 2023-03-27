#!/usr/bin/env python3

import numpy as np
from pylab import *
import math
import os
import subprocess
import time
import shutil

originalFile='originalFiles/spotweldModelOriginal.k'

class CoorDesc:
    
    paramCurrent=[] # of order of dim:  paramCurrent = [x1,x2,...,xn]
    paramDelta=[]   # of order of dim:  paramDelta = [dx1,dx2,...,dxn] 
    paramRange=[]   #=[ [x1Min,x1Max], [x2Min,x2Max], ...[xnMin,xnMax] ]
    
    valueCurrent=None # Store the value of function v=f(paramCurrent)
    
    
    indx=0                  # indx in range  0,1,...,(dim-1)
    direction=0             # 0:Neutral, +1:Positive, -1:Negative
    growthStatus=1          # 0:Neutral, +1:Expanding, -1:Contracting
    currentError=None       # Current Value of Error
    nextRunParamList=None   # List of param for next Run: [ [run1_x1,run1_x2], [run2_x1,run2_x2] ]
    
    winner=None
    caseID=None
    
    maxItrinSameIndx=7
    ItrinSameIndx=0
    
    def __init__(self,paramCurrent,paramRange,paramDelta=[],nCPUs=2):
        initialDeltaRatio=20   # if not input by user: dx1=(x1Max-x1Min)/20, dx2=(x2Max-x2Min)/20, ...
        self.dim = len(paramCurrent)
        self.nCPUs = nCPUs
        self.paramCurrent=paramCurrent
        self.paramRange=paramRange
        
        if len(paramDelta)==0:
            self.paramDelta=[ (x[1]-x[0])/initialDeltaRatio for x in paramRange]
        #-------End of Initialization-------------------------------------------
        
    
    def changeIndx(self):
        if self.indx<len(self.paramCurrent)-1:
            self.indx +=1
        else:
            self.indx=0
        
        self.direction=0
        self.growthStatus=1
        self.ItrinSameIndx=0
    
    def updateByError(self,errorList,valueList=[]):
        if len(valueList)==0:
            valueList=errorList
        
        if not self.currentError:   # 1st iteration
            self.currentError = errorList[0]
            self.valueCurrent=valueList[0]
            return
        
        if self.currentError <= min(errorList):
            winner='P0'
        elif errorList[0] <= errorList[1]:
            winner='P1'
            self.valueCurrent=valueList[0]
        else:
            winner='P2'
            self.valueCurrent=valueList[1]
        
        if self.direction==0 and self.growthStatus==1:      #case1
            caseID='case1'
            if winner=='P0':
                self.paramDelta[self.indx] /=2
                self.growthStatus=-1
            elif winner=='P1':
                self.paramCurrent=self.nextRunParamList[0].copy()
                self.currentError=errorList[0]
                self.direction=-1
                self.paramDelta[self.indx] *=2   
            elif winner=='P2':
                self.paramCurrent=self.nextRunParamList[1].copy()
                self.currentError=errorList[1]
                self.direction=1
                self.paramDelta[self.indx] *=2 
            
        elif self.direction==1 and self.growthStatus==1:    #case2
            caseID='case2'
            if winner=='P0':
                self.paramDelta[self.indx] /=2
                self.changeIndx()
            elif winner=='P1':
                self.paramCurrent=self.nextRunParamList[0].copy()
                self.currentError=errorList[0]
                self.changeIndx()
                
            elif winner=='P2':
                self.paramCurrent=self.nextRunParamList[1].copy()
                self.currentError=errorList[1]
                self.paramDelta[self.indx] *=2
        
        elif self.direction==-1 and self.growthStatus==1:    #case3
            caseID='case3'
            if winner=='P0':
                self.paramDelta[self.indx] /=2
                self.changeIndx()
            elif winner=='P1':
                self.paramCurrent=self.nextRunParamList[0].copy()
                self.currentError=errorList[0]
                self.changeIndx()
            elif winner=='P2':
                self.paramCurrent=self.nextRunParamList[1].copy()
                self.currentError=errorList[1]
                self.paramDelta[self.indx] *=2
        elif self.direction==0 and self.growthStatus==-1:    #case4
            caseID='case4'
            if winner=='P0':
                self.paramDelta[self.indx] /=2
            elif winner=='P1':
                self.paramCurrent=self.nextRunParamList[0].copy()
                self.currentError=errorList[0]
                self.changeIndx() 
            elif winner=='P2':
                self.paramCurrent=self.nextRunParamList[1].copy()
                self.currentError=errorList[1]
                self.changeIndx()
                
        elif self.direction==1 and self.growthStatus==-1:    #case5
            caseID='case5'
            if winner=='P0':
                self.paramDelta[self.indx] /=2
            elif winner=='P1':
                self.paramCurrent=self.nextRunParamList[0].copy()
                self.currentError=errorList[0]
                self.changeIndx()
            elif winner=='P2':
                self.paramCurrent=self.nextRunParamList[1].copy() 
                self.currentError=errorList[1]
                self.paramDelta[self.indx] /=2
                self.changeIndx()                
        elif self.direction==-1 and self.growthStatus==-1:    #case6
            caseID='case6'
            if winner=='P0':
                self.paramDelta[self.indx] /=2
            elif winner=='P1':
                self.paramCurrent=self.nextRunParamList[0].copy()
                self.currentError=errorList[0]
                self.changeIndx() 
            elif winner=='P2':
                self.paramCurrent=self.nextRunParamList[1].copy()
                self.currentError=errorList[1]
                self.paramDelta[self.indx] /=2
                self.changeIndx()
        
        self.caseID=caseID
        self.winner=winner
        
        self.ItrinSameIndx +=1
        if self.ItrinSameIndx >self.maxItrinSameIndx:
            self.changeIndx()
        
        return 

    def NextRunList(self):       
        runParamList=[]
        
        if not self.currentError:   # 1st iteration
            runParamList.append(self.paramCurrent)
            self.nextRunParamList=runParamList
            return runParamList
        #------------------------------------------

        for iRun in range(self.nCPUs):
            runNew = self.paramCurrent.copy()
            
            if self.direction==0:
                dx = (-1)**(iRun+1) * self.paramDelta[self.indx]
            else:
                dx = (self.direction)* self.paramDelta[self.indx]** \
                (self.growthStatus*(iRun+1)) 

            runNew[self.indx] += dx
            
            runNew[self.indx] = max(runNew[self.indx],
                                  self.paramRange[self.indx][0])
                                  
            runNew[self.indx] = min(runNew[self.indx],
                                  self.paramRange[self.indx][1])
            runParamList.append(runNew)
            
            
        print('runParamList:',runParamList)
        
        self.nextRunParamList=runParamList
        return runParamList        

def giveFullFolderAddress(subFolder,geoConfig):
    return 'runFolder/'+geoConfig+'/'+subFolder
    
def read_DispLoad_Data(filePath,isExperimental=False):
    nLineStart = 0 if (isExperimental) else 2
    loadScale = 0.001 if (isExperimental) else 1
    delim = '\t' if (isExperimental) else ','
    
    f = open(filePath, 'r+')
    Lines = f.readlines()
    f.close()
    
    disp=[]
    load=[]
    
    iLine=0
    for line in Lines:
        if iLine>= nLineStart:
            splited=line.split(delim)
            disp.append(float(splited[0]))
            load.append(loadScale*float(splited[1]))
        iLine+=1
    return disp,load

def generateSubFolderList(nCores):
    # Generate Folder List:----------------------
    subFolderList=[]

    for iSubFolder in range(1,nCores+1):
        subFolderList.append('Processor_'+str(iSubFolder)+"/")
    #--------------------------------------------
   
    return subFolderList

def materialCarts(varVecIn):
    
    matCoeffs=[0]*11
    cartJumpIndx=[8]    # Means that cart2 starts from indx #8.     Example: [3,5]:cart1:[0,1,2], cart2[3,4], cart3[5:end]
    cartLines=[7,9]     # Two carts in line #7 & line #9
    coeffVarIndx=[2,9]
    
    
    matCoeffs[0]=40200143    #[0]    mid
    matCoeffs[1]=1.0850E-6   #[1]    ro
    matCoeffs[2]=0.9         #[2]    e***          
    matCoeffs[3]=0.15        #[3]    pr
    matCoeffs[4]=0.3         #[4]    sigy   
    matCoeffs[5]=0           #[5]    eh       
    matCoeffs[6]=0           #[6]    dt
    matCoeffs[7]=0           #[7]    tfail

    matCoeffs[8]=0           #[8]    efail
    matCoeffs[9]=0.05        #[9]    sigLcAx***                      
    matCoeffs[10]=0          #[10]   sigLcAtu
    
    if (not len(varVecIn)==len(coeffVarIndx) or not (len(varVecIn)==len(cartLines)) ):
        pritn('Error! Number of variable coefficients is inconsistent!')
        exit()
    
    
    carts=[]
    indx=0
    indxVecIn=0
    
    for i in range(len(cartJumpIndx)+1):
        if i<len(cartJumpIndx):
            iMax=cartJumpIndx[i]
        else:
            iMax=len(matCoeffs)
        
        cart=[]
        while indx < iMax:
            if indx in coeffVarIndx:
                cart.append(varVecIn[indxVecIn])
                indxVecIn+=1
            else:
                cart.append(matCoeffs[indx])
            indx+=1
        
        carts.append(cart)
        
    return carts,cartLines
    
def checkIfRunFinished(folder):
    filePath = folder+'pyStatus.txt'
    f = open(filePath, 'r+')
    Lines = f.readlines()
    f.close()
    if len(Lines)==0:
        return False
    else:
        if (Lines[0]=='Running'):
            return False
        else:
            return True

def runFolders(subFolderList,nRuns,geoConfigs):
    currentFolder=os.getcwd()+"/"

    nFolders=nRuns
    isFinished=[False]*nFolders
    
    
    
    
    for config in geoConfigs:
        folderList=[giveFullFolderAddress(subFolder,config) for subFolder in subFolderList]
        
        for iFolder in range(1,nFolders+1):
            targetFolder=folderList[iFolder-1]
            os.chdir(currentFolder+targetFolder)
            subprocess.Popen(["python.exe", "pyRunDyna.py"],creationflags=subprocess.CREATE_NEW_CONSOLE)
            os.chdir(currentFolder)

    print('Running...')
    time.sleep(1)

    for config in geoConfigs:
        folderList=[giveFullFolderAddress(subFolder,config) for subFolder in subFolderList]
        for iRun in range(1,nFolders+1):
            targetFolder=folderList[iFolder-1]
            while not isFinished[iRun-1]:
                isFinished[iRun-1] = checkIfRunFinished(targetFolder)
    
    time.sleep(2)
    
    return 0

def modifyMaterialFile(origin,desti,cartsTxt,cartLines):
    f = open(origin, 'r+')
    Lines = f.readlines()
    f.close()
    
    f = open(desti, 'w')
    iLine=0
    indxCart=0
    
    for line in Lines:
        
        iLine +=1
        
        if ( indxCart>=len(cartsTxt) ) or not (iLine==cartLines[indxCart]):
            f.write(line)
        else:
            f.write(cartsTxt[indxCart])
            indxCart+=1
            
    f.close()
    
    return 0

def createCartsTxt(carts):
    
    cartsTxt=[]
    for iCart in range(len(carts)):
        vecIn=carts[iCart]
        n=len(vecIn)
        txtCart=""
        for i in range(n):
            L=len(str(vecIn[i]))
            if L<10:
                for j in range(10-L):
                    txtCart +=" "
                txtCart += str(vecIn[i])
            else:
                txtCart += '{: .3E}'.format(vecIn[i])
            
        txtCart +="\n"
        cartsTxt.append(txtCart)
    return cartsTxt

def assignParameterValues(paramRunList,subFolderList,geoConfigs):
    
    for config in geoConfigs:
        for iFolder in range(len(paramRunList)):
            destinationFile=giveFullFolderAddress(subFolderList[iFolder],config)+"spotweldModel.k"
            carts,cartLines=materialCarts(paramRunList[iFolder])
            cartsTxt=createCartsTxt(carts)
            modifyMaterialFile(originalFile,destinationFile,cartsTxt,cartLines)
    return 0

def readSimulationResults(subFolderList,nRuns,geoConfigs):
    '''
    peakListSim=[   [
                        [dispPeak_LSP,loadPeak_LSP],[dispPeak_LSV,loadPeak_LSV],[dispPeak_LSV,loadPeak_CPH] \core1
                    ],
                    [
                        [dispPeak_LSP,loadPeak_LSP],[dispPeak_LSV,loadPeak_LSV],[dispPeak_LSV,loadPeak_CPH] \core2
                    ],
                    [
                        [dispPeak_LSP,loadPeak_LSP],[dispPeak_LSV,loadPeak_LSV],[dispPeak_LSV,loadPeak_CPH] \core3
                    ],
                    ...
                ]
    '''
    peakListSim=[]
    
    
    for iFolder in range(nRuns):
        peakCurrentFolder=[]
        for config in geoConfigs:
            folder=giveFullFolderAddress(subFolderList[iFolder],config)
            filePath = folder + 'LOAD_vs_Disp'
            disp,load = read_DispLoad_Data(filePath)

            dispPeak=disp[load.index(max(load))]
            loadPeak=max(load)
            peakCurrentFolder.append([dispPeak,loadPeak])
            
        peakListSim.append(peakCurrentFolder)        
    return peakListSim

def calcErrors(givenXY_List,targetXYs):

    errors=[]
    for iError in range(len(givenXY_List)):
        PeaksCurrentRun=givenXY_List[iError]
        err=0
        for iConfig in range(len(PeaksCurrentRun)):
            err += ((PeaksCurrentRun[iConfig][0]-targetXYs[iConfig][0])/targetXYs[iConfig][0])**2 + \
               ((PeaksCurrentRun[iConfig][1]-targetXYs[iConfig][1])/targetXYs[iConfig][1])**2
        errors.append(err/len(PeaksCurrentRun))
    return errors

def calcNextCoeffs(coeffValueOld,coeffDelta,errors,learningRate):
    nPar=len(coeffValueOld)
    
    derivatives=[None]*nPar
    coeffValueNew=[None]*nPar
    
    for iPar in range(nPar):
        derivatives[iPar]=(errors[2*iPar+2]-errors[2*iPar+1])/(2*coeffDelta[iPar])
        coeffValueNew[iPar]=coeffValueOld[iPar] - learningRate*coeffDelta[iPar]*derivatives[iPar]
    
    
    return coeffValueNew

#def recordResults(itr,coeffValue,peaks,error):
def recordResults(itr,OptObj,geoConfigs,paramRunList,peakListExp,peakListSim):
    outputFile='Results/numerical/convergenceHistory'
    
    nConfigs=len(geoConfigs)
    
    coeffs=OptObj.paramCurrent
    peaks=OptObj.valueCurrent
    error=OptObj.currentError
    
    nPar=len(coeffs)
    
    print('-----------------------------------------------')
    print('\n iteration:{}\n'.format(itr+1))
    print('Best Coefficients:',coeffs)
    print('Best Curve Peaks:',peaks)
    print('Lowest Error:',error)
    print('-------------------')
    print('Parameters Set:',paramRunList)
    print('Experimental Peaks:',peakListExp)
    print('Simulation Peaks:',peakListSim)
    
    # Writing into recording file:
    if itr==0:
        writtenLine='Iteration,'
        
        for iPar in range(nPar):
            writtenLine+='parameter #{},'.format(iPar+1)
            writtenLine+='TotalError,'
        

        for iConfig in range(nConfigs):
            config=geoConfigs[iConfig]
            
            writtenLine+='{}_dispPeak,'.format(config)
            writtenLine+='{}_loadPeak'.format(config) 
            writtenLine+="\r" if iConfig==len(geoConfigs)-1 else ","  
        
        
        f = open(outputFile, 'w')
        f.write(writtenLine)
        f.write('-'*len(writtenLine))
        f.write('\r')
        f.close()
    
    writtenLine='{},'.format(itr+1)

    for iPar in range(nPar):
        writtenLine+=('{:.5E},'.format(coeffs[iPar]) )
    
    writtenLine+=('{:.5E},'.format(error) )
    
    for iConfig in range(nConfigs):
    
            writtenLine+='{:.5E},'.format(peaks[iConfig][0])
            writtenLine+='{:.5E}'.format(peaks[iConfig][1])
            writtenLine+="\r" if iConfig==len(geoConfigs)-1 else ","
    
    f = open(outputFile, 'a')
    f.write(writtenLine)
    f.close()

    return 0

def readExpResults(expResFolder,geoConfigs):
    '''
    peakExpList=[
                    [dispPeak_LSP,loadPeak_LSP], 
                    [dispPeak_LSV,loadPeak_LSV],
                    [dispPeak_LSV,loadPeak_CPH], 
                    ...
                ]
    '''

    fileName='LOAD_vs_Disp'
    peakExpList=[]
    for config in geoConfigs:
        fileFullName=expResFolder+config+"/"+fileName
        disp,load=read_DispLoad_Data(fileFullName,isExperimental=True)
        dispPeak=disp[load.index(max(load))]
        loadPeak=max(load)
        peakExpList.append([dispPeak,loadPeak])
    return peakExpList


nCores=2 
geoConfigs=['LSP','LSV']
expResFolder='Results/experimental/'

peakListExp = readExpResults(expResFolder,geoConfigs)

coeffValue=[0.1,0.03] # Initial Guess
coeffLim1=[0.01,2]
coeffLim2=[0.005,0.7]

OptObj=CoorDesc(coeffValue,[coeffLim1,coeffLim2])
errorList=[]


subFolderList = generateSubFolderList(nCores)
itrMax=int(input('Enter Number of Iterations:'))
errorOK=1e-5


for itr in range(0,itrMax):

    paramRunList=OptObj.NextRunList()
    nRuns=len(paramRunList)
    
    assignParameterValues(paramRunList,subFolderList,geoConfigs)
    runFolders(subFolderList,nRuns,geoConfigs)
    peakListSim=readSimulationResults(subFolderList,nRuns,geoConfigs)
    errorList=calcErrors(peakListSim,peakListExp)   
    OptObj.updateByError(errorList,peakListSim)
    recordResults(itr,OptObj,geoConfigs,paramRunList,peakListExp,peakListSim)

print('-----------------------------------------------')
print('Finished!')


 

