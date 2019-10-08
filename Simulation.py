# -*- coding: utf-8 -*-
"""
Created on Wed April 24 08:18:01 2019
Lasted modified on August 31, 2019

authors: Yihong Mauro, John Mauro, Collin Wilkinson

"""

import Constants 
import Temperature
import numpy as np
import copy

class MyException(Exception):
    pass


class Simulation:
      
    def __init__(self, inParameterFileName):
        inParameters = {}
        inFileName = inParameterFileName
        #paramFileName= input("Please enter the file name for input parameters:  ") # Se_test5       
        try:
            fhParam = open(inFileName, 'r')
        except IOError:
            print('Cannot open the input parameter file: ', inFileName)
        else:
            line_param = fhParam.readline()            
            while line_param:       
                if line_param.startswith("#"): # ignore the lines that start with the comment sign #
                    continue
                if (';' in line_param):
                    line_param = line_param.replace(";","") # remove ";"
                elif ('#' in line_param):    
                    line_param = line_param[0:line_param.find('#')] # remove everything after the comment sign #
                 
                line_param_splitted = line_param.split() #splits on all whitespace if left unspecified
                inParameters["%s"%(line_param_splitted[0])] = line_param_splitted[2]
                line_param = fhParam.readline()        
        fhParam.close()  
 
        self.mParameters = inParameters
        print("mParameters = ")
        print(self.mParameters)

        requiredParams = ["SimulationName", "InitialTemperature", "SimulationTime", 
                         "OutputResolution", "NumBasins", "NumAtoms", "AtomicMass", 
                         "ObservationTime", "HopDistance", "PrintPhaseDistributions", 
                         "TemperaturePath"]
        requiredPositiveParams = ["InitialTemperature", "OutputResolution", "NumBasins", 
                                 "NumAtoms", "AtomicMass", "HopDistance"]
        
        for par in requiredParams:
            if par not in self.mParameters:
                raise MyException("%s is required in your input file.\n"%(par))
        
        for parp in requiredPositiveParams:
            if float(self.mParameters[parp]) <= 0:
                raise MyException("The value of %s needs to be positive.\n"%(parp))
   
        self.mSimulationName = self.mParameters["SimulationName"]       
        self.mInitialTemperature = float( self.mParameters["InitialTemperature"] )        
        self.mSimulationTime = float( self.mParameters["SimulationTime"] )
        self.mOutputResolution = float( self.mParameters["OutputResolution"] )
        self.mNumBasins = int(float( self.mParameters["NumBasins"]))
        self.mNumAtoms = int(float( self.mParameters["NumAtoms"] ))
        self.mAtomicMass = float( self.mParameters["AtomicMass"] ) / Constants.kAvogadrosNumber / 1000.0
        self.mObservationTime = float( self.mParameters["ObservationTime"] )
        self.mHopDistance = float( self.mParameters["HopDistance"] )
        self.mPrintPhaseDistributions = float( self.mParameters["PrintPhaseDistributions"] )


        ###################################################
        #    Temperature path     
        ###################################################
  
        if ( self.mParameters["TemperaturePath"] ==  "Fixed" ):
                      
            if "Temperature" not in self.mParameters:
                raise MyException("Temperature is required in the input parameter file for Fixed TemperaturePath. \n")
            elif self.mParameters["Temperature"] <= 0:
                raise MyException("Temperature is required to be positive. \n")
            else:
                self.mTemperaturePath = Temperature.FixedTemperature( float( self.mParameters["Temperature"] ) )
        
        elif ( self.mParameters["TemperaturePath"] == "Linear" ):
            requiredparms = ["FinalTemperature", "TotalCoolingTime"]   
            for i in range(len(requiredparms)):        
                if requiredparms[i] not in self.mParameters:
                    raise MyException("%s is required in your input file. \n"%(requiredparms[i]))
    
            for i in range(len(requiredparms)):
                if float(self.mParameters[requiredparms[i]]) <= 0:
                    raise MyException("The value of %s needs to be positive. \n"%(requiredparms[i]))
     
            self.mTemperaturePath = Temperature.LinearTemperature( self.mInitialTemperature, 
                                                  float( self.mParameters["FinalTemperature"] ), 
                                                  float( self.mParameters["TotalCoolingTime"] ) )
               
        elif ( self.mParameters["TemperaturePath"] == "FromFile" ):
            if FileName not in self.mParameters:
                raise MyException("FileName is required in your input file. \n")
                
            self.mTemperaturePath = Temperature.FromFileTemperature( self.mInitialTemperature, self.mParameters["FileName"] )
   
        elif ( self.mParameters["TemperaturePath"] == "QuenchHoldLinear" ):
            requiredparms = ["QuenchTemperature", "HoldTime", "FinalTemperature", "CoolingTime"]
    
            for i in range(len(requiredparms)):        
                if requiredparms[i] not in self.mParameters:
                    raise MyException("%s is required in your input file. \n"%(requiredparms[i]))
    
            for i in range(len(requiredparms)):
                if float(self.mParameters[requiredparms[i]]) <= 0:
                    raise MyException("The value of %s needs to be positive. \n"%(requiredparms[i]))
    
            self.mTemperaturePath = Temperature.QuenchHoldLinearTemperature(float( self.mParameters["QuenchTemperature"] ),
                                                                float( self.mParameters["HoldTime"] ),
                                                                float( self.mParameters["FinalTemperature"] ),
                                                                float( self.mParameters["CoolingTime"] ) )
        
        elif ( self.mParameters["TemperaturePath"] == "Exponential" ):
            requiredparms = ["FinalTemperature", "DecayTime"]
       
            for i in range(len(requiredparms)):        
                if requiredparms[i] not in self.mParameters:
                    raise MyException("%s is required in your input file. \n"%(requiredparms[i]))
    
            for i in range(len(requiredparms)):
                if float(self.mParameters[requiredparms[i]]) <= 0:
                    raise MyException("The value of %s needs to be positive. \n"%(requiredparms[i]))
    
            self.mTemperaturePath = Temperature.ExponentialTemperature(self.mInitialTemperature,
                                                      float( self.mParameters["FinalTemperature"] ),
                                                      float( self.mParameters["DecayTime"] ) )
   
        elif ( self.mParameters["TemperaturePath"] == "Zigzag" ): 
            requiredparms = ["LowTemperature", "HighTemperature", "FinalTemperature", 
                             "InitialCoolingTime", "ReheatingTime", "FinalCoolingTime",
                             "LowHoldTime", "HighHoldTime"]
            requiredpositive = ["LowTemperature", "HighTemperature", "FinalTemperature"]
    
            for i in range(len(requiredparms)):        
                if requiredparms[i] not in self.mParameters:
                    raise MyException("%s is required in your input file. \n"%(requiredparms[i]))
    
            for i in range(len(requiredpositive)):
                if float(self.mParameters[requiredpositive[i]]) <= 0:
                    raise MyException("The value of %s needs to be positive. \n"%(requiredpositive[i]))
    
            self.mTemperaturePath = Temperature.ZigzagTemperature(self.mInitialTemperature,
                                                 float( self.mParameters["LowTemperature"] ),
                                                 float( self.mParameters["HighTemperature"] ),
                                                 float( self.mParameters["FinalTemperature"] ),
                                                 float( self.mParameters["InitialCoolingTime"] ),
                                                 float( self.mParameters["ReheatingTime"] ),
                                                 float( self.mParameters["FinalCoolingTime"] ),
                                                 float( self.mParameters["LowHoldTime"] ),
                                                 float( self.mParameters["HighHoldTime"] ) )
        else:
            raise MyException ( "The specified TemperaturePath is not implemented." )


        #############################################################################
        #     Read in the energy landscape files: 
        #          - take the file names from the input parameter file 
        #          - five file names ended with .min, .prp, .deg, .trs, and .eig   
        #############################################################################

        requiredFilenames = ["MinimumFile", "TransitionFile", "PropertyFile", "EigenvalueFile", "DegeneracyFile"]
        for file in requiredFilenames:   
            if file not in self.mParameters:
                raise MyException("%s is required in your input file. \n"%(file))

        #mMinimumNumbers = np.zeros(mNumBasins, dtype = int)
        self.mMinimumNumbers = np.arange(self.mNumBasins, dtype = int)
        self.mEnergyMatrix = 1000.0*np.ones((self.mNumBasins, self.mNumBasins), dtype = float)
        self.mEffectiveEnergyMatrix = np.zeros((self.mNumBasins, self.mNumBasins), dtype = float)
        self.mPropertyVector = np.zeros(self.mNumBasins)
        self.mFrequencyMatrix =  np.zeros((self.mNumBasins, self.mNumBasins), dtype = float)
        self.mLogDegeneracy = np.zeros(self.mNumBasins, dtype = float)
        
        self.mNumMetabasins = 1
        self.mMetabasinRates = list(np.zeros(self.mNumMetabasins))
        self.mMetabasinRates[0] = np.zeros((self.mNumBasins, self.mNumBasins), dtype = float) # np array (various sizes) in each list element
        self.mMetabasinRank = np.zeros(self.mNumMetabasins, dtype = int)
        self.mMetabasinRank[0] = self.mNumBasins
        self.mMetabasinIndices = list(np.zeros(self.mNumMetabasins, dtype = int))
        self.mMetabasinIndices[0] = np.arange(self.mNumBasins, dtype = int) # np array (various sizes) in each list element
        
        #self.mOldMetabasinProbabilities = np.zeros( self.mNumMetabasins )
        self.mMetabasinProbabilities = np.zeros( self.mNumMetabasins)       
        
       
        ##### read in Minimum file (.min file): [basin#, minimum_energy]
        try:
            fh_min = open(self.mParameters["MinimumFile"], "r")
        except IOError:
            print('cannot open', fh_min)
        else:
            lines_min = fh_min.readlines()
            if len(lines_min) != self.mNumBasins:
                raise MyException("Length of minimum file does not match number of basins.")
            
            for i in range(self.mNumBasins):
                line_min_splitted = (lines_min[i]).split()
                
                if len(line_min_splitted) != 2 :
                    raise MyException("Minimum file should have two columns." )
    
                self.mMinimumNumbers[i] = int(float(line_min_splitted[0]))
                self.mEnergyMatrix[i, i] = float(line_min_splitted[1]) * Constants.kElectronVolt#/Constants.kBoltzmannsConstant
                
            self.mBiasEnergy = (-1.0) *np.min(np.diagonal(self.mEnergyMatrix))

            row,col = np.diag_indices(self.mEnergyMatrix.shape[0])
            self.mEnergyMatrix[row,col] = np.diagonal(self.mEnergyMatrix) + self.mBiasEnergy
            #mEnergyMatrix = mEnergyMatrix + mBiasEnergy
          
        fh_min.close()


        ##### read in Property file (.prp file): [basin#, molar_volume]
        try:
            fh_prp = open(self.mParameters["PropertyFile"], "r")
        except IOError:
            print('cannot open', fh_prp)
        else:
            lines_prp = fh_prp.readlines()
            #print(lines_prp)
            if len(lines_prp) != self.mNumBasins:
                raise MyException("Length of property file does not match number of basins.")
    
            for i in range(self.mNumBasins):
                line_prp_splitted = (lines_prp[i]).split()
                
                if len(line_prp_splitted) != 2 :
                    raise MyException("Property file should have two columns." )
                           
                theIndex = np.where(self.mMinimumNumbers == int(float(line_prp_splitted[0])))[0]
                self.mPropertyVector[theIndex] = float( line_prp_splitted[1] )
           
        fh_prp.close()
    
 
        ##### read in Degeneracy file (.deg file): [basin#, log(gi)]
        try:
            fh_deg = open(self.mParameters["DegeneracyFile"], "r")
        except IOError:
            print('cannot open', fh_deg)
        else:
            lines_deg = fh_deg.readlines()
            if len(lines_deg) != self.mNumBasins:
                raise MyException("Length of degeneracy file does not match number of basins.")
    
            for i in range(self.mNumBasins):
                line_deg_splitted = (lines_deg[i]).split()
               
                if len(line_deg_splitted) != 2 :
                    raise MyException("Degeneracy file should have two columns." )
                           
                theIndex = np.where(self.mMinimumNumbers == int(float(line_deg_splitted[0])))[0]
                self.mLogDegeneracy[theIndex] = float( line_deg_splitted[1] )
            
    
        fh_deg.close()
    
    
        ##### read in Transition file (.trs file): 
        #[transition point#, energy at the transition point, basins that are connected the the transition point (left to right)]
        try:
            fh_trs = open(self.mParameters["TransitionFile"], "r")
        except IOError:
            print('cannot open', fh_trs)
        else:
            line_trs = fh_trs.readlines()
            for i in range(len(line_trs)):
                line_trs_splitted = (line_trs[i]).split()

                if len(line_trs_splitted) != 4 :
                    raise MyException("Transition file should have four columns." )
                           
                theFirstIndex = np.where(self.mMinimumNumbers == int(float(line_trs_splitted[2])))[0]
                theSecondIndex = np.where(self.mMinimumNumbers == int(float(line_trs_splitted[3])))[0]
    
                if ( theSecondIndex < theFirstIndex ):
                    theTempIndex = theFirstIndex
                    theFirstIndex = theSecondIndex
                    theSecondIndex = theTempIndex
    
                theTransitionEnergy = float( line_trs_splitted[1] ) * Constants.kElectronVolt + self.mBiasEnergy   
                self.mEnergyMatrix[theFirstIndex, theSecondIndex] = theTransitionEnergy
 
        fh_trs.close()

    
        ##### read in Eigenvalue file (.eig file): [basin#1, basin#2, eigenvalue]
        try:
            fh_eig = open(self.mParameters["EigenvalueFile"], "r")
        except IOError:
            print('cannot open', fh_eig)
        else:
            line_eig = fh_eig.readlines()
            #print(line_eig)
            for i in range(len(line_eig)):
                line_eig_splitted = (line_eig[i]).split()      
                #print("line_eig_splitted = ")
                #print(line_eig_splitted)
                if (len(line_eig_splitted) != 3):
                    raise MyException("Eigenvalue file should have three columns." )
    
                theFirstIndex = np.where(self.mMinimumNumbers == int(float(line_eig_splitted[0])))[0]
                theSecondIndex = np.where(self.mMinimumNumbers == int(float(line_eig_splitted[1])))[0]
                theEigenvalue = float( line_eig_splitted[2] ) 
                self.mFrequencyMatrix[theFirstIndex, theSecondIndex] = ( theEigenvalue / self.mAtomicMass )**0.5 / ( 2.0 * np.pi )
                print(self.mFrequencyMatrix)
        
        fh_eig.close()


        #############################################################################
        #    More initialization
        #############################################################################
        #a = np.array([[1, 1, 0, 0],
        #  [0, 0, 1, 1],
        #  [0, 0, 0, 0]])
        #idx = np.asarray(np.where(a == 1))
        #a[[*idx]]= -1000
        #    indices = np.asarray(np.where(mEnergyMatrix == 1000))
        #    mEnergyMatrix[[*indices]] = -1000000000

        
        ### fill in more transition energies
        for i in range(self.mNumBasins):
            for j in range(i+1, self.mNumBasins):
                if ( self.mEnergyMatrix[i,j] == 1000 ):
                    self.mEnergyMatrix[i,j] = -1000000000
                    
                    for k in range(self.mNumBasins):
                        if ( k != i ): 
                            if ( ( self.mEnergyMatrix[i,k] != 1000 ) and ( self.mEnergyMatrix[i,k] > self.mEnergyMatrix[i,j] ) ):
                                self.mEnergyMatrix[i,j] = self.mEnergyMatrix[i,k]
        
                        if ( k != j ):
                            if ( ( self.mEnergyMatrix[k,j] != 1000 ) and ( self.mEnergyMatrix[k,j] > self.mEnergyMatrix[i,j] ) ):
                                 self.mEnergyMatrix[i,j] = self.mEnergyMatrix[k,j]
          
                    if ( self.mEnergyMatrix[i,j] == -1000000000 ): 
                        raise MyException( "Could not find transition!\n" )

        # create a symmetric 2D array using the upper triangle
        i_lower = np.tril_indices(self.mNumBasins, -1)
        self.mEnergyMatrix[i_lower] = self.mEnergyMatrix.T[i_lower]  
     
        for i in range(self.mNumBasins):
            for j in range(self.mNumBasins):
                if ( ( i != j ) and self.mFrequencyMatrix[i, j] == 0 ):
                    for k in range(self.mNumBasins):
                        if ( (self.mFrequencyMatrix[i, k] != 0) and (self.mFrequencyMatrix[k, j] != 0) ):
                            self.mFrequencyMatrix[i, j] = self.mFrequencyMatrix[i, k]
     
        for i in range(self.mNumBasins): 
            self.mEffectiveEnergyMatrix[i,i] = self.mEnergyMatrix[i,i] - Constants.kBoltzmannsConstant * self.mInitialTemperature * self.mLogDegeneracy[i]
            for j in range(i + 1, self.mNumBasins):
                self.mEffectiveEnergyMatrix[i,j] = self.mEnergyMatrix[i,j] - Constants.kBoltzmannsConstant * self.mInitialTemperature * self.mLogDegeneracy[j]    
                self.mEffectiveEnergyMatrix[j,i] = self.mEnergyMatrix[j,i] - Constants.kBoltzmannsConstant * self.mInitialTemperature * self.mLogDegeneracy[i]

        self.mPhaseSpaceVector = np.exp(-np.diagonal(self.mEffectiveEnergyMatrix) / ( Constants.kBoltzmannsConstant * self.mInitialTemperature ) )
        thePartitionFunction = np.sum(self.mPhaseSpaceVector)
        self.mPhaseSpaceVector /= thePartitionFunction + 1.0e-999
               
        print("\n__init__() is done!\n")

   
    def EquilibriumMetabasins(self, inTemperature):
        
        for i in range( self.mNumMetabasins ):
            thePartitionFunction = 0            

            for j in range( self.mMetabasinRank[i] ):
                theIndex = self.mMetabasinIndices[i][j]
                self.mPhaseSpaceVector[theIndex] = np.exp( -self.mEffectiveEnergyMatrix[theIndex][theIndex] / ( Constants.kBoltzmannsConstant * inTemperature ) )
                thePartitionFunction += self.mPhaseSpaceVector[theIndex]
                
            for j in range( self.mMetabasinRank[i] ):
                theIndex = self.mMetabasinIndices[i][j]
                self.mPhaseSpaceVector[theIndex] /= thePartitionFunction + 1.0e-999
                self.mPhaseSpaceVector[theIndex] *= self.mMetabasinProbabilities[i]

    
    def MetabasinProbability(self, inIndex):
        outProbability = 0
        if ( inIndex < self.mNumMetabasins ):
            for i in range(self.mMetabasinRank[inIndex]):
                outProbability += self.mPhaseSpaceVector[self.mMetabasinIndices[inIndex][i]]
        return outProbability


    def CalculateInterMetabasinRates(self, inTemperature):              
        for i in range( self.mNumMetabasins ):
            for j in range( self.mNumMetabasins ):
                self.mInterMetabasinRates[i][j] = 0
          
        for i in range(self.mNumMetabasins):
            if ( self.mMetabasinProbabilities[i] != 0 ):
                for j in range(self.mNumMetabasins):
                    if ( self.mUseMetabasinTransition[i][j] == True ):
                        if ( i != j ):
                            for k in range(self.mMetabasinRank[i]):
                                for l in range(self.mMetabasinRank[j]):
                                    theRow = self.mMetabasinIndices[i][k]
                                    theCol = self.mMetabasinIndices[j][l]
                                    theRate = self.mPhaseSpaceVector[theRow] * self.mFrequencyMatrix[theRow, theCol] * np.exp ( -( self.mEffectiveEnergyMatrix[theRow, theCol] - self.mEnergyMatrix[theRow, theRow] ) / ( Constants.kBoltzmannsConstant * inTemperature ) ) / self.mMetabasinProbabilities[i]
                                    self.mInterMetabasinRates[i][j] += theRate
                                    self.mInterMetabasinRates[i][i] -= theRate


    def DetermineMetabasins(self, inTemperature, inTimeStep):
        theThreshold = 1.0 / inTimeStep

        theCount = 0
        theRateMatrix = np.zeros( (self.mNumBasins, self.mNumBasins), dtype = float )
        theMetabasinMatrix = np.zeros( (self.mNumBasins, self.mNumBasins), dtype = bool )
        theUseMetabasinArray = np.zeros( self.mNumBasins, dtype = bool )

        for i in range(self.mNumBasins):
            for j in range(self.mNumBasins):
                if ( i != j ):
                    	theRateMatrix[i, j] = self.mFrequencyMatrix[i,j] * np.exp(-( self.mEffectiveEnergyMatrix[i,j] - self.mEnergyMatrix[i,i] ) / ( Constants.kBoltzmannsConstant * inTemperature ) )
      
        for i in range(self.mNumBasins):
            for j in range(i + 1, self.mNumBasins):
                if ( ( theRateMatrix[i, j] > theThreshold ) and ( theRateMatrix[j, i] > theThreshold ) ):
                    	theMetabasinMatrix[i, j] = True
                else:
                    	theMetabasinMatrix[i, j] = False

        while True:                     
            for i in range(self.mNumBasins):
                for j in range(i + 1, self.mNumBasins):
                    if ( theMetabasinMatrix[i, j] == True):
                        for k in range(i + 1, j):
                            if ( theMetabasinMatrix[k, j] == True):
                                theMetabasinMatrix[i, k] = True
                        
            theUseMetabasinArray = np.ones( self.mNumBasins, dtype = bool ) # all elements are True           
            theCount = 0
            for i in range(self.mNumBasins):
                if ( theUseMetabasinArray[i] == True):
                    theCount += 1
                    for j in range(i + 1, self.mNumBasins):
                        if ( theMetabasinMatrix[i, j] == True ): 
                            theCount += 1
                            if ( theUseMetabasinArray[j] == True):
                                theUseMetabasinArray[j] = False
                            
            if ( theCount == self.mNumBasins ):
                break
  
        while True:
            theMetabasinsFinished = True      
            theUseMetabasinArray = np.ones( self.mNumBasins, dtype = bool ) # all elements are True

            self.mNumMetabasins = 0
            theCount = 0
            
            for i in range(self.mNumBasins):
                if ( theUseMetabasinArray[i] == True ):
                    self.mNumMetabasins += 1
                    theCount += 1
                                         
                    for j in range (i + 1, self.mNumBasins):
                        if ( theMetabasinMatrix[i, j] == True ):
                            theCount += 1
                            
                            if ( theUseMetabasinArray[j] == True ):
                                theUseMetabasinArray[j] = False
              
            self.mMetabasinRates = list(np.zeros( self.mNumMetabasins, dtype = float ))
            self.mMetabasinRank = np.zeros( self.mNumMetabasins, dtype = int )
            self.mMetabasinIndices = list(np.zeros( self.mNumMetabasins, dtype = int ))
            
            theMetabasinIndex = 0
            theBasinIndex = 0           
            for i in range(self.mNumBasins):
                if ( theUseMetabasinArray[i] == True ):
                    self.mMetabasinRank[theMetabasinIndex] = 1
                    
                    for j in range(i + 1, self.mNumBasins):
                        if ( theMetabasinMatrix[i, j] == True ):
                            self.mMetabasinRank[theMetabasinIndex] += 1
	                           
                    self.mMetabasinRates[theMetabasinIndex] = np.zeros( (self.mMetabasinRank[theMetabasinIndex], self.mMetabasinRank[theMetabasinIndex]), dtype = float )
                    self.mMetabasinIndices[theMetabasinIndex] = np.zeros ( self.mMetabasinRank[theMetabasinIndex], dtype = int )           
                    self.mMetabasinIndices[theMetabasinIndex][0] = i
                    
                    theBasinIndex = 1
                    for j in range(i + 1, self.mNumBasins):
                        if ( theMetabasinMatrix[i, j] == True ):
                            self.mMetabasinIndices[theMetabasinIndex][theBasinIndex] = j
                            theBasinIndex += 1
                            
                    theMetabasinIndex += 1
             
            theRow = 0
            theCol = 0
            for i in range( self.mNumMetabasins ):
                for j in range( self.mMetabasinRank[i] ):
                    for k in range( self.mMetabasinRank[i] ):
                        if ( j != k ):
                            theRow = self.mMetabasinIndices[i][j]
                            theCol = self.mMetabasinIndices[i][k]
                            self.mMetabasinRates[i][j][k] = self.mFrequencyMatrix[theRow][theCol] * np.exp ( -( self.mEffectiveEnergyMatrix[theRow][theCol] - self.mEnergyMatrix[theRow][theRow] ) / ( Constants.kBoltzmannsConstant * inTemperature ) )
          
            self.mMetabasinProbabilities = np.zeros( self.mNumMetabasins )
            self.mOldMetabasinProbabilities = np.zeros( self.mNumMetabasins )
            
            for i in range(self.mNumMetabasins):
                self.mMetabasinProbabilities[i] = self.MetabasinProbability( i )
                self.mOldMetabasinProbabilities[i] = self.mMetabasinProbabilities[i]
    

            self.mUseMetabasinTransition = np.zeros( (self.mNumMetabasins, self.mNumMetabasins), dtype = bool )
            self.mInterMetabasinRates = np.zeros( (self.mNumMetabasins, self.mNumMetabasins ), dtype = float)
                         
            for i in range( self.mNumMetabasins ):
                for j in range( self.mNumMetabasins ):
                    if ( i == j ):
                        self.mUseMetabasinTransition[i, j] = False
                    else:
                        self.mUseMetabasinTransition[i, j] = True
    
            self.CalculateInterMetabasinRates( inTemperature )

  
            for i in range(self.mNumMetabasins):
                for j in range(i + 1,  self.mNumMetabasins):
                    if theMetabasinsFinished == True:
                        if ( ( self.mInterMetabasinRates[i][j] > theThreshold ) and ( self.mInterMetabasinRates[j][i] > theThreshold ) ):
                            theMetabasinsFinished = False
                            theMetabasinMatrix[self.mMetabasinIndices[i][0]][self.mMetabasinIndices[j][0]] = True
     
            if ( theMetabasinsFinished == True):
                break

        theTotalRank = np.sum(self.mMetabasinRank)
        assert(theTotalRank == self.mNumBasins)
        
        return theThreshold


    def MasterEquationSolver(self):
                
        theFileName= input("Please enter the name of the output file:  ") 
        #theFileName = self.mSimulationName + "_pt01.txt" #77989s
        theFile= open(theFileName , 'w')

        theAverageProperty = 0
        theEquilibriumProperty = 0
        theEquilibriumEntropy = 0
        theOldTemperature = 0
        theOldEnergy = 0
        thePartitionFunction = 0
        theEnergy = (-1.) * self.mBiasEnergy
        theAverageLogDegeneracy = 0
        theEquilibriumLogDegeneracy = 0
        theEquilibriumEnergy = (-1.) * self.mBiasEnergy
        thePhaseTrace = 0
        theEquilibriumPhaseTrace = 0
        theTemperature = self.mTemperaturePath.Evaluate(0)
        theTemperatureStep = 0
        theDerivative = 0
        theTime = 0
        theTimeStep = 0
        theLastOutput = 0
        theThreshold = 0

        
        self.mEquilibriumPhase = copy.deepcopy(self.mPhaseSpaceVector)
        theEnergy += np.sum(self.mPhaseSpaceVector * np.diagonal(self.mEnergyMatrix))

        theEnergy *= Constants.kAvogadrosNumber
        
        theOldTemperature = self.mInitialTemperature
        theOldEnergy = theEnergy
        
        outstring = "#Time_s\tTemp_K\tProp\tEq_Prop\tU_kJ/mol\t"
        outstring += "EqU_kJ/mol\tlng\teq_lng\t"
        outstring += "Cp_J/mol-K\tNum_MB\tTrPhase\tTrPhEq\t1/T"

        theFile.write(outstring)
        
        if ( self.mPrintPhaseDistributions ):
            for i in range(self.mNumBasins):
                theFile.write("\t%d"%(i + 1))
                    
            for i in range(self.mNumBasins):
                theFile.write("\t%d"%(i + 1))
    
        theFile.write("\n")
        
        #theTimeStep = 0.01
        
        print("starting MasterEquationSolver\n")


        while ( theTime < self.mSimulationTime ):
            theTimeStep = self.mTemperaturePath.SafeTimeStep( theTime )
            #theTimeStep = 0.01
            if ( theTimeStep == 0 ):
                theTimeStep = 1.0e-4 * self.mSimulationTime
                
            if ( theTime + theTimeStep > self.mSimulationTime ):
                theTimeStep = self.mSimulationTime - theTime
                                                
#            theTemperature = self.mTemperaturePath.Evaluate( theTime )
#            #print("theTime = %8.5f\n"%(theTime))
#            theTemperature = 0
#            theCount = 0
#            for  t  in np.arange(theTime, theTime + theTimeStep, 0.025 * theTimeStep ):
#                theTemperature += self.mTemperaturePath.Evaluate ( t )
#                theCount += 1
#                
#            theTemperature /= theCount
#            theTime += theTimeStep

            theTime += theTimeStep

            #if ( theTime > self.mSimulationTime ):
            #    theTime = self.mSimulationTime

            theTemperature = self.mTemperaturePath.Evaluate ( theTime )       
           
            for i in range(self.mNumBasins):
                self.mEffectiveEnergyMatrix[i,i] = self.mEnergyMatrix[i,i] - Constants.kBoltzmannsConstant * theTemperature * self.mLogDegeneracy[i]
                
                for j in range(i + 1, self.mNumBasins):
                    self.mEffectiveEnergyMatrix[i,j] = self.mEnergyMatrix[i,j] - Constants.kBoltzmannsConstant * theTemperature * self.mLogDegeneracy[j]
                    self.mEffectiveEnergyMatrix[j,i] = self.mEnergyMatrix[j,i] - Constants.kBoltzmannsConstant * theTemperature * self.mLogDegeneracy[i]
            
            theThreshold = self.DetermineMetabasins( theTemperature, theTimeStep )
            self.EquilibriumMetabasins ( theTemperature )

            #print("theTime = %8.5f, theTemperature = %8.5f, mNumMetabasins = %d\n"%(theTime, theTemperature, self.mNumMetabasins))          
                       
#            self.mUseMetabasinTransition = np.zeros( ( self.mNumMetabasins, self.mNumMetabasins ), dtype = bool) # all True
#            self.mInterMetabasinRates = np.zeros( (self.mNumMetabasins, self.mNumMetabasins ), dtype = float)
#
#            for i in range(self.mNumMetabasins):
#                for j in range(self.mNumMetabasins):
#                    if ( i == j ):
#                        self.mUseMetabasinTransition[i][j] = False
#                    else:
#                        self.mUseMetabasinTransition[i][j] = True
#                        
#            self.CalculateInterMetabasinRates( theTemperature )
            
            
            if (self.mNumMetabasins > 1):
                theDeltaProb = 0
                
                for i in range(self.mNumMetabasins):
                    for j in range(i + 1, self.mNumMetabasins):
                        if ( self.mUseMetabasinTransition[i, j] ):
                            if ( ( self.mInterMetabasinRates[i][j] > theThreshold ) and ( self.mInterMetabasinRates[j,i] < theThreshold ) ):
                                theDeltaProb = ( self.mMetabasinProbabilities[i] * self.mInterMetabasinRates[i][j] - self.mMetabasinProbabilities[j] * self.mInterMetabasinRates[j][i] ) / ( self.mInterMetabasinRates[i][j] + self.mInterMetabasinRates[j][i] )
                                self.mMetabasinProbabilities[i] -= theDeltaProb
                                self.mMetabasinProbabilities[j] += theDeltaProb
                                
                                self.EquilibriumMetabasins ( theTemperature )
                                
                                for k in range(self.mNumMetabasins):
                                    self.mUseMetabasinTransition[i,k] = False
                                    self.mUseMetabasinTransition[k,i] = False
                                    
                                self.CalculateInterMetabasinRates ( theTemperature )
                                
                            elif ( ( self.mInterMetabasinRates[j][i] > theThreshold ) and ( self.mInterMetabasinRates[i][j] < theThreshold ) ):
                                theDeltaProb = ( self.mMetabasinProbabilities[j] * self.mInterMetabasinRates[j][i] - self.mMetabasinProbabilities[i] * self.mInterMetabasinRates[i][j] ) / ( self.mInterMetabasinRates[j][i] + self.mInterMetabasinRates[i][j] )
                                self.mMetabasinProbabilities[j] -= theDeltaProb
                                self.mMetabasinProbabilities[i] += theDeltaProb
                                
                                self.EquilibriumMetabasins( theTemperature )
                                
                                for k in range(self.mNumMetabasins):
                                    self.mUseMetabasinTransition[j][k] = False
                                    self.mUseMetabasinTransition[k][j] = False
                                    
                                self.CalculateInterMetabasinRates( theTemperature )
                                
                            elif ( ( self.mInterMetabasinRates[j][i] > theThreshold ) and ( self.mInterMetabasinRates[i][j] > theThreshold ) ):
                                print("THIS SHOULD NEVER HAPPEN!\n")
                                
                theEulerTimeStep = theTimeStep * 0.01
                theEulerTime = 0
                
                while ( theEulerTime < theTimeStep ):
                    theEulerTime += theEulerTimeStep
                    self.mOldMetabasinProbabilities = copy.deepcopy(self.mMetabasinProbabilities)
                    
                    #for i in range(self.mNumMetabasins):
                    #    self.mOldMetabasinProbabilities[i] = self.mMetabasinProbabilities[i]
                        
                    for i in range(self.mNumMetabasins):
                        for j in range(self.mNumMetabasins):
                            if ( self.mUseMetabasinTransition[i][j] and ( i != j ) ):
                                theDerivative = self.mInterMetabasinRates[j][i] * self.mOldMetabasinProbabilities[j] - self.mInterMetabasinRates[i][j] * self.mOldMetabasinProbabilities[i]
                                self.mMetabasinProbabilities[i] += theDerivative * theEulerTimeStep
                                
                    self.EquilibriumMetabasins( theTemperature )
                    self.CalculateInterMetabasinRates( theTemperature )
                    
            if ( ( theTime - theLastOutput ) >= self.mOutputResolution ):
                theLastOutput = theTime
                thePhaseTrace = 0
                theEquilibriumPhaseTrace = 0
                theAverageProperty = 0
                theEnergy = (-1)*self.mBiasEnergy
                theEquilibriumEnergy = (-1)*self.mBiasEnergy
                theAverageLogDegeneracy = 0
                theEquilibriumLogDegeneracy = 0
                                 
                for i in range(self.mNumBasins):
                    theAverageProperty += self.mPhaseSpaceVector[i] * self.mPropertyVector[i]
                    theEnergy += self.mPhaseSpaceVector[i] * self.mEnergyMatrix[i][i]
                    
                    theInnerSum = 0
                    for j in range(self.mNumBasins):
                        if ( j != i ):
                            theInnerSum += self.mPhaseSpaceVector[j] * np.exp( self.mLogDegeneracy[i] )
                            
                    theAverageLogDegeneracy += self.mPhaseSpaceVector[i] * theInnerSum
                    
                theAverageLogDegeneracy = np.log( theAverageLogDegeneracy )
                thePhaseTrace = np.sum(self.mPhaseSpaceVector)  #self.PhaseTrace()
                #print("Time = %ds, #MBs = %d, PhaseTrace = %8.5f, %5.2f, Completed\n"%(theTime, self.mNumMetabasins, thePhaseTrace, theTime / self.mSimulationTime * 100.0))
                #assert ( thePhaseTrace <= 1.01 ) and ( thePhaseTrace >= 0.99 ) 
                
                thePartitionFunction = 0
                for i in range(self.mNumBasins):
                    self.mEquilibriumPhase[i] = np.exp( self.mLogDegeneracy[i] - self.mEnergyMatrix[i][i] / ( Constants.kBoltzmannsConstant * theTemperature ) )
                    thePartitionFunction += self.mEquilibriumPhase[i]
                    
                #for i in range(self.mNumBasins):
                #    self.mEquilibriumPhase[i] /= thePartitionFunction
                
                self.mEquilibriumPhase /= thePartitionFunction + 1.0e-999
                    
                theEquilibriumPhaseTrace = np.sum(self.mEquilibriumPhase) #self.EquilibriumPhaseTrace()
                
                theEquilibriumProperty = 0
                theEnsembleEntropy = 0
                theEquilibriumEntropy = 0
                theEquilibriumLogDegeneracy = 0
                
                for i in range(self.mNumBasins):
                    theEquilibriumProperty += self.mEquilibriumPhase[i] * self.mPropertyVector[i]
                    theEquilibriumEnergy += self.mEquilibriumPhase[i] * self.mEnergyMatrix[i][i]
                    theInnerSum = 0
                    
                    for j in range(self.mNumBasins):
                        if ( j != i ):
                            theInnerSum += self.mEquilibriumPhase[j] * np.exp( self.mLogDegeneracy[j] )
                            
                    theEquilibriumLogDegeneracy += self.mEquilibriumPhase[i] * theInnerSum
                    
                    if ( self.mPhaseSpaceVector[i] > 0 ):
                        theEnsembleEntropy -= self.mPhaseSpaceVector[i] * np.log( self.mPhaseSpaceVector[i] )
                        
                    if ( self.mEquilibriumPhase[i] > 0 ):
                         theEquilibriumEntropy -= self.mEquilibriumPhase[i] * np.log( self.mEquilibriumPhase[i] )
                         
                theEquilibriumLogDegeneracy = np.log( theEquilibriumLogDegeneracy )
                theEnergy *= Constants.kAvogadrosNumber
                theEquilibriumEnergy *= Constants.kAvogadrosNumber
                                             
                theFile.write("%8.5f\t%8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\t%d\t%8.5e\t%8.5e\t%8.5e"
                              %(theTime, theTemperature, theAverageProperty, theEquilibriumProperty, 
                                theEnergy, theEquilibriumEnergy, theAverageLogDegeneracy, 
                                theEquilibriumLogDegeneracy, 
                                ( theEnergy - theOldEnergy ) / ( theTemperature - theOldTemperature ), 
                                self.mNumMetabasins, thePhaseTrace, theEquilibriumPhaseTrace, 1.0 / theTemperature))

                print("theTime = %8.5f, theTemperature = %8.5f, mNumMetabasins = %d"%(theTime, theTemperature, self.mNumMetabasins))          
                print("theAverageProperty = %8.5e, theEnergy = %8.5e"%(theAverageProperty, theEnergy))
                print("mEquilibriumPhase = ")
                print(self.mEquilibriumPhase)
                print("\n")
                
                if ( self.mPrintPhaseDistributions ):
                    for i in range(self.mNumBasins):
                        theFile.write( "\t%8.5e"%(self.mPhaseSpaceVector[i]))
                       
                    for i in range(self.mNumBasins):
                        theFile.write("\t%8.5e"%(self.mEquilibriumPhase[i]))
                       
                theFile.write("\n")
    
                theOldTemperature = theTemperature
                theOldEnergy = theEnergy

        print("Final phase space distribution vector:\n")
        print(self.mPhaseSpaceVector)
        print("\n\n")
        
        theFile.close()
 
                     