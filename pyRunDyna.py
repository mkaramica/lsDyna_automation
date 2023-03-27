#!/usr/bin/env python3

import numpy as np
from pylab import *
import os

def readSimulationResults(filePath):

    f = open(filePath, 'r+')
    Lines = f.readlines()
    f.close()
    
    strCheckLoad="force resultants"
    strCheckTime="output at time ="
    
    disp=[]
    load=[]
    iLine=0
    for line in Lines:
        iLine +=1
        if strCheckLoad in line:
            newLine=line[36:-1]
            splited=newLine.split("   ")
            load.append(-float(splited[2]))
        if strCheckTime in line:
            newLine=line[18:-1]
            t=float(newLine)
            disp.append(0.005*t)
    #caseID=input('Enter CaseID:')       
    plt.plot(disp,load)
    plt.xlabel('Disp (mm)')
    plt.ylabel('Load (kN)')
    #plt.title('caseID:'+caseID)
    plt.grid(True)
    plt.show()
    
    return disp,load


def deletePreviousFiles(exceptions):
    allFiles = [f for f in os.listdir('.') if os.path.isfile(f)]
    for f in allFiles:
        if not f in exceptions:
            os.remove(f)
    return 0

def recordStatus(status):
    fileName='pyStatus.txt'
    f = open(fileName, 'w')
    f.write(status)
    f.close()
    return 0

def saveResults(outputFile):
    disp,load=readSimulationResults('spcforc')
    
    '''
    #Remove negative oscillations:
    indxPeak=load.index(max(load))
    indx1stNegative=-1
    for iPoint in range(len(disp)):
        if iPoint>indxPeak and load[iPoint]<0:
            disp=disp[0:iPoint]
            load=load[0:iPoint]
            break
    
    disp.append(disp[-1])
    load.append(0)
    '''
    f = open(outputFile, 'w')
    f.write('disp,Load\n' )
    f.write('--------------------------\n' )
    for i in range(len(disp)):
        f.write('%8.6e,%8.6e\n' % (disp[i],load[i]) )
    f.close()
    
    return 0

outputFile='LOAD_vs_Disp'

mainfileList = ['Control_Inc.k',
              'MID48260201_TPO.k',
              'spotweldModel.k',
              'mainKey.key',
              'runcommand.bat',
              'pyStatus.txt',
              os.path.basename(__file__)]   # Current Python Script


deletePreviousFiles(mainfileList)
recordStatus('Running')
os.system('runcommand.bat')

if os.path.exists("spcforc"):
    recordStatus('Finished')
    saveResults(outputFile)
    mainfileList.append(outputFile)
else:
    recordStatus('Failed')


#deletePreviousFiles(mainfileList)
    
