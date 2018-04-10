import os
import glob
from xml.dom import minidom
import Averager
import EdfFile as edf
import CenterOfRotation
import numpy as np


class FastTomoExperiment:
    def __init__(self, xmlFilename, output):
        self.dirname = os.path.dirname(xmlFilename)
        self.projectionName = os.path.dirname(xmlFilename)
        self.outputFolder = output
        self.xmlfile = xmlFilename
        self.scanRange = 0
        self.numberOfProjections = 0
        self.wfEveryFloor = True
        self.wfInterval = True
        self.wfIntervalWithinAFrame = 0
        self.ref_on = 0
        self.numberOfDarkFields = 0
        self.numberOfProjections = 0
        self.projectionNumber180Degree = 0
        self.numberFlatField = 0
        self.pixelSize = 46
        self.halfAcquisition = False
        self.halfAcquiSide = 'left'
        self.ringFilterActivated = 0
        self.paganinMode = 0
        self.paganinLength = 1000
        self.startLine = 1
        self.energy=0
        self.beamline='ID17'
        self.nameExp='XXXX'
        self.machineCurrentStart=0
        self.machineCurrentStop=0
        self.xmldoc = minidom.parse(self.xmlfile)
        self.defineCorrectValues()

        self.averageDone = False
        self.normalizedTomography=False
        self.cor = -1

    def setHalfAcquisition(self, side):
        self.halfAcquisition = 1
        self.halfAcquiSide = side

    def setPaganin(self, length):
        self.paganinMode = 1
        self.paganinLength = length

    def setRingFilter(self):
        self.ringFilterActivated = 1

    def setLineToReconstruct(self, beg, end):
        self.startLine = beg
        self.endLine = end

    def setMissingPoints(self, ddict):

        if (ddict["ringFilter"] == True):
            self.setRingFilter()
        if (ddict["paganin"] == True):
            self.setPaganin(ddict["paganinLength"])

        if (ddict["wholeStack"] == False):
            print
            self.setLineToReconstruct(int(ddict["startLine"]), int(ddict["endLine"]))

        if (ddict["halfAcquisition"] == True):
            self.setHalfAcquisition(ddict["halfAcquiSide"])


    def defineCorrectValues(self):

        for node in self.xmldoc.getElementsByTagName("scanRange"):
            self.scanRange = float(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("machineCurrentStart"):
            self.machineCurrentStart= float(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("machineCurrentStop"):
            self.machineCurrentStop= float(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("beamline"):
            self.beamline= str(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("nameExp"):
            self.nameExp= str(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("disk"):
            self.pathExp= str(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("tomo_N"):
            self.numberOfProjections = int(node.childNodes[0].toxml())

            if self.scanRange == 180:
                self.projectionNumber180Degree = self.numberOfProjections

            else:
                self.projectionNumber180Degree = self.numberOfProjections / 2.

        for node in self.xmldoc.getElementsByTagName("ref_On"):
            self.ref_on = int(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("dark_N"):
            self.numberOfDarkFields = int(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("ref_N"):
            self.numberFlatField = int(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("DIM_1"):
            self.widthProjection = int(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("DIM_2"):
            self.heightProjection = int(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("scanName"):
            self.scanName = str(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("pixelSize"):
            self.pixelSize = float(node.childNodes[0].toxml())

        for node in self.xmldoc.getElementsByTagName("scanName"):
            self.projectionName = (node.childNodes[0].toxml())


    def createAverageWfandDf(self):
        refBeg = glob.glob(self.dirname + '/ref*_0000.edf')
        if len(refBeg)>0:
            self.outputFileFFNameBeg = self.dirname + '/refForHST0000.edf'
            vBeg = Averager.Averager(refBeg, self.outputFileFFNameBeg, option=0)

        if(self.ref_on<self.numberOfProjections):
            cpt=self.ref_on
            while cpt<self.numberOfProjections:
                textref = '%4.4d' % cpt
                refBeg = glob.glob(self.dirname + '/ref*_' + textref + '.edf')

                if len(refBeg>0):
                    self.outputFileNameFFEnd = self.dirname + '/refForHST' + textref + '.edf'
                    vBeg = Averager.Averager(refBeg, self.outputFileNameFFEnd, option=1)
                cpt+=self.ref_on
        else:
            textref = '%4.4d' % self.ref_on
            refBeg = glob.glob(self.dirname + '/ref*_' + textref + '.edf')
            if len(refBeg)>0:
                self.outputFileNameFFEnd = self.dirname + '/refForHST' + textref + '.edf'
                vBeg = Averager.Averager(refBeg, self.outputFileNameFFEnd, option=0)

        darkFile = self.dirname + '/darkend0000.edf'
        data = edf.EdfFile(darkFile).GetData(0)
        data /= self.numberOfDarkFields
        self.darkOutputFile = self.dirname + '/darkForHST0000.edf'

        filetoWrite = edf.EdfFile(self.darkOutputFile, access='wb+')
        filetoWrite.WriteImage({}, data, Append=0, DataType='FloatValue')

        self.averageDone = True


    def findCenterOfRotation(self):
        self.cor = -1
        if (self.averageDone == False) :
            self.createAverageWfandDf()


        self.project0Name = self.dirname + '/' + self.scanName + '0000.edf'
        text180 = '%4.4d' % self.projectionNumber180Degree
        self.project180Name = self.dirname + '/' + self.scanName + text180 + '.edf'
        edfFile0 = edf.EdfFile(self.project0Name, access='rb')
        project0 = edfFile0.GetData(0)

        edfFile0 = edf.EdfFile(self.project180Name, access='rb')
        project180 = edfFile0.GetData(0)

        edfFile0 = edf.EdfFile(self.outputFileFFNameBeg, access='rb')
        ff = edfFile0.GetData(0)

        meanEnergy = np.mean(ff, 1)
        lineofMaxExnergy = np.argmax(meanEnergy)
        ff[ff == 0] = 1
        project0[project0 == 0] = 0.0000000001
        project180[project180 == 0] = 0.0000000001
        project0 = (np.true_divide(project0, ff))
        project180 = (np.true_divide(project180, ff))
        project0[np.isnan(project0) == True] = 0.0000000001
        project180[np.isnan(project180) == True] = 0.0000000001

        project0 = -np.log(project0)
        project180 = -np.log(project180)

        if (self.halfAcquisition):
            #we need to cut the projection in two
            # but we need to know the side
            halfWidth = self.widthProjection / 2.

            if (self.halfAcquiSide == 'left'):
                project180 = project180[:, 0:int(halfWidth)]
                project180 = np.fliplr(project180)
                project0 = project0[:, 0:int(halfWidth)]

            else:
                project180 = project180[:, int(halfWidth):int(halfWidth) * 2]
                project180 = np.fliplr(project180)
                project0 = project0[:, int(halfWidth):int(halfWidth) * 2]

        else:
            project180 = np.fliplr(project180)


        numberOfLinesForAverage = int(self.heightProjection / 4.)
        centralLine = int(self.heightProjection / 2.)
        if (self.scanRange ==180):
            self.cor = CenterOfRotation.findCorLudo(project0, project180, centralLine, numberOfLinesForAverage)
        else:
            print '360 degree acquisition multiple cor calculations'
            self.cor = CenterOfRotation.findCorLudo(project0, project180, centralLine, numberOfLinesForAverage)


        if (self.halfAcquiSide == 'right'):
            self.cor += int(halfWidth)

        return self.cor

    def askMissingQuestions(self):
        self.halfAcquisition = 0
        condition = False
        while (condition == False):
            var = raw_input("Acquisition made in Half Tomography? (y,n) : ")
            print
            'You entered var :' + var
            if (var == "Y") or (var == "y") or (var == "N") or (var == "n"):
                condition = True

        if ((var == 'Y') or (var == 'y')):
            self.halfAcquisition = 1
            condition = False
            while (condition == False):
                varCor = raw_input("Center of rotation on the left ? (y,n) : ")
                if ((varCor == 'Y') or (varCor == 'y') or (varCor == 'N') or (varCor == 'n')):
                    condition = True

            if ((varCor == 'Y') or (varCor == 'y')):
                self.halfAcquiSide = 'left'
            else:
                self.halfAcquiSide = 'right'

        else:
            self.halfAcquisition = 0


    def printOut(self):
        print 'scanRangeFound ' + str(self.scanRange)
        print 'number of projection ' + str(self.numberOfProjections)
        print 'projection 180  ' + str(self.projectionNumber180Degree)
        print 'angular step is therefore ' + str(self.scanRange / self.numberOfProjections)
        print 'Number of Dark Field : ' + str(self.numberOfDarkFields)
        print 'Number of white Field : ' + str(self.numberFlatField)
        print 'Width Of The Projections : ' + str(self.widthProjection)
        print 'Height Of The Projections : ' + str(self.heightProjection)
        print ' ref_on : ' + str(self.ref_on)



if __name__ == "__main__":
    framesDirectory = '/Volumes/ID17/map3/Papyrus/Blue_Banana_real/Blue_banana_real_70keV_037_'
    scanParameterFile = glob.glob(framesDirectory + '/*.xml')
    scanParameterFile.sort()
    print scanParameterFile

    tomoExp = FastTomoExperiment(scanParameterFile[0],'/')
    tomoExp.printOut()
    corFound=tomoExp.findCenterOfRotation()
    print 'Cor Found at ' +str(corFound)
# print acquisitionNode


