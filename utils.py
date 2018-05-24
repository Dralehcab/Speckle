from random import randint
import glob
import os
import spytIO
import numpy as np




def generateArrayOfRandomNumbers(nbPoints,deb,end):
    listOfPoint=[]
    if(deb<end):
        while(len(listOfPoint)<nbPoints):
            number=randint(deb,end)
            if not(number in listOfPoint):
                listOfPoint.append(number)

    return listOfPoint


def spytMkDir(folder):
    if not (os.path.exists(folder)):
        os.mkdir(folder)


def convert32BitFolder(inputFolder,outputFolder):
    listofFiles=glob.glob(inputFolder+'/*.edf')
    for filename in listofFiles:
        base=os.path.basename(filename)
        data=spytIO.openImage(filename)
        dataToStore=np.asarray(data,np.float32)
        outputfilename=outputFolder+'/'+base
        print ('Copying \n'+str(filename)+' \n in \n'+str(outputfilename))
        spytIO.saveEdf(dataToStore,outputfilename)


def convert32BitopticalFlowFolder(inputFolder):
    print('Convert entire Experiment')

    dxInput=inputFolder+'/dx/'
    dxOutput=inputFolder+'/dx32/'
    spytMkDir(dxOutput)
    convert32BitFolder(dxInput,dxOutput)

    dyInput = inputFolder + '/dy/'
    dyOutput = inputFolder + '/dy32/'
    spytMkDir(dyOutput)
    convert32BitFolder(dyInput, dyOutput)

    phiInput = inputFolder + '/phi/'
    phiOutput = inputFolder + '/phi32/'
    spytMkDir(phiOutput)
    convert32BitFolder(phiInput, phiOutput)

    phi2Input = inputFolder + '/phi2/'
    phi2Output = inputFolder + '/phi2_32/'
    spytMkDir(phi2Output)
    convert32BitFolder(phi2Input, phi2Output)

    phi3Input = inputFolder + '/phi2/'
    phi3Output = inputFolder + '/phi2_32/'
    spytMkDir(phi3Output)
    convert32BitFolder(phi3Input, phi3Output)




if __name__ == "__main__":
    list=generateArrayOfNumbers(10,0,58)
    print(list)

    inputFolder='/Volumes/ID17/speckle/md1097/id17/SpeckleSacroIliaque/OpticalFlowFloor0/'
    convert32BitopticalFlowFolder(inputFolder)