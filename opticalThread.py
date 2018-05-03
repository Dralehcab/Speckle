from threading import Thread
import spytIO
import corrections as corr
import OpticalFlow
import os



def createFolder(folder):
    if not (os.path.exists(folder)):
        os.mkdir(folder)


class OpticalFlowSolver(Thread):

    def __init__(self, listOfDictionaries,listOfProjections,outputFolder):
        Thread.__init__(self)
        self.listOfDictionnaries = listOfDictionaries
        self.listOfProjection=listOfProjections
        self.output=outputFolder
        self.createFolders()

    def createFolders(self):
        self.dxFolder = self.output + '/dx/'
        self.dyFolder = self.output + '/dy/'
        self.phiFolder = self.output + '/phi/'
        self.phi2Folder = self.output + '/phi2/'
        self.phi3Folder = self.output + '/phi3/'

        createFolder(self.dxFolder)
        createFolder(self.dyFolder)
        createFolder(self.phiFolder)
        createFolder(self.phi2Folder)
        createFolder(self.phi3Folder)

    def run(self):
        for numeroProjection in self.listOfProjection:
            numeroProjection=int(numeroProjection)
            print('Processing '+str(numeroProjection))
            projectionFiles = []
            referencesFiles = []
            darkFieldFiles = []
            for dict in self.listOfDictionnaries:
                projectionFiles.append(dict['projections'][numeroProjection])
                referencesFiles.append(dict['references'][0])
                darkFieldFiles.append(dict['darkField'])


            print(projectionFiles)

            Is = spytIO.openSeq(projectionFiles)
            Ir = spytIO.openSeq(referencesFiles)
            df = spytIO.openSeq(darkFieldFiles)
            #Is, Ir = corr.registerImagesBetweenThemselves(Is, Ir)
            result = OpticalFlow.processProjectionSet(Is, Ir, df)

            dx = result['dx']
            dy = result['dy']
            phi = result['phi']
            phi2 = result['phi2']
            phi3 = result['phi3']

            textProj='%4.4d'%numeroProjection
            spytIO.saveEdf(dx.real, self.dxFolder + '/dx' + textProj + '.edf')
            spytIO.saveEdf(dy.real, self.dyFolder + '/dy' + textProj + '.edf')
            spytIO.saveEdf(phi.real, self.phiFolder + '/phi' + textProj + '.edf')
            spytIO.saveEdf(phi2.real, self.phi2Folder + '/phi2' + textProj + '.edf')
            spytIO.saveEdf(phi3.real, self.phi3Folder + '/phi3' + textProj + '.edf')






