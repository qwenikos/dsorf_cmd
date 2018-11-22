import os
import sys
sys.path.append('./lib')
from CP_features_extaction_module import *
from nMersDB import *
from fn_module import *


baseOutputDir="defaultOutput/"
currentPath=os.path.dirname(os.path.realpath(sys.argv[0]))+"/"

###DefaultValues Start
outputDir=""
configFileName="conf/sORFconfig.cfg"
numOfProcess=1 ##for parallel processing of module
mode=1
bypassSignalPeptide="N"
signalPeptideBypaseMultiplier=0 ##select 0 to not bypassing signlal peptide zone
simulateLength=""
###DefaultValues End

if len(sys.argv)==9:
    #print "argument are ok "
    inputData   =sys.argv[1]
    outputDir       =sys.argv[2]
    numOfProcess    =sys.argv[3]
    mode            =sys.argv[4]
    startingPos     =sys.argv[5]
    signalPeptideBypaseMultiplier=sys.argv[6]
    configFileName=sys.argv[7]
    input_type      =sys.argv[8] ##0 for seq 1 for filename
    
if len(sys.argv)==10:
    #print "argument are ok "
    inputData   =sys.argv[1]
    outputDir       =sys.argv[2]
    numOfProcess    =sys.argv[3]
    mode            =sys.argv[4]
    startingPos     =sys.argv[5]
    signalPeptideBypaseMultiplier=sys.argv[6]
    configFileName=sys.argv[7]
    simulateLength=sys.argv[8]
    input_type      =sys.argv[9] ##0 for seq 1 for file name
    
if input_type=="0":
    tempFileName="tempInputFile.fa"
    tempFile=open(tempFileName,"w")
    tempFile.write(">tempSorfId\n")
    tempFile.write(inputData+"\n")
    tempFile.close()
    inputData="tempInputFile.fa"
    

if int(startingPos)<7 and (not (mode==2)):
    print "cannot continue with this mode due to upstream region length. Select  mode 2 (Coding Combisition only)"
    exit()

if bypassSignalPeptide=="N":
    signalPeptideBypaseMultiplier=0
elif bypassSignalPeptide=="Y":
    signalPeptideBypaseMultiplier=3

print signalPeptideBypaseMultiplier


###read config file params
confDict=readConfigFile(configFileName)
outputStartingDir=confDict["outputStartingDir"]
outputStartingDir=currentPath+outputStartingDir
minLenLim=confDict["minLenLim"]
#########################
if (not outputStartingDir=="") and  (not outputStartingDir[-1:]=="/"):
    outputStartingDir=outputStartingDir+"/"
    
if (not outputDir=="") and (not outputDir[-1:]=="/"):
    outputDir=outputDir+"/"

workingDir=outputStartingDir+outputDir
print workingDir

command="mkdir -p "+workingDir
os.system(command)

#########################

command="python src/call_predict.py "+inputData+" "+workingDir+" "+str(signalPeptideBypaseMultiplier)+" "+str(startingPos)+" "+str(mode)+" "+str(numOfProcess)+" "+minLenLim+" "+configFileName +" "+inputData+" "+ simulateLength
print command
os.system(command)
    
print "END"



