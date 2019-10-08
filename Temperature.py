# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 07:54:09 2019
Last modified: Apr 21, 2019

@author: Yihong Mauro and Yihong Mauory

"""


########################################################################################
#     Abstract Temperature class: an arbitrary temperature path in time, T(t)
########################################################################################

from abc import ABC, abstractmethod
import numpy as np

class Temperature(ABC):    
    @abstractmethod
    def Evaluate(self):
        pass

    @abstractmethod
    def SafeTimeStep(self):
        pass    


########################################################################################
#     LinearTemperature class: inherite from Temperature class
########################################################################################

class LinearTemperature (Temperature):
    def __init__(self, inInitialTemperature, inFinalTemperature, inTotalCoolingTime):
        self.mInitialTemperature = inInitialTemperature
        self.mFinalTemperature = inFinalTemperature
        self.mTotalCoolingTime = inTotalCoolingTime        
          
    def Evaluate(self, inTime): 
        if (inTime > self.mTotalCoolingTime ):
            return self.mFinalTemperature        
        if ( inTime < 0 ):
            return self.mInitialTemperature        
        return self.mInitialTemperature + ( self.mFinalTemperature - self.mInitialTemperature ) / self.mTotalCoolingTime * inTime
    
    def SafeTimeStep(self, inTime, inTemperatureStep = 0.01):
        if ( inTime >= self.mTotalCoolingTime ):
            return 0
        if ( inTime <= 0 ):
            return 0;
        return inTemperatureStep*self.mTotalCoolingTime/abs(self.mFinalTemperature-self.mInitialTemperature)


########################################################################################
#     ExponentialTemperature class: inherits from Temperature class
########################################################################################

class ExponentialTemperature(Temperature):
    def __init__(self, inInitialTemperature, inFinalTemperature, inDecayTime):
        self.mInitialTemperature = inInitialTemperature
        self.mFinalTemperature = inFinalTemperature
        self.mDecayTime = inDecayTime        
         
    def Evaluate(self, inTime): 
        if ( inTime < 0 ):
            return self.mInitialTemperature
        return self.mInitialTemperature-(self.mFinalTemperature-self.mInitialTemperature)*(1.0-np.exp(-inTime/self.mDecayTime))

    def SafeTimeStep(self, inTime, inTemperatureStep = 0.001 ):
        if ( inTime <= 0 ):
            return 0
        return inTemperatureStep*1.0e-12/abs(self.Evaluate(inTime)-self.Evaluate(inTime+1.0e-12))
       
        
########################################################################################
#     QuenchHoldLinearTemperature class: inherits from Temperature class
########################################################################################

class QuenchHoldLinearTemperature(Temperature):
    def __init__(self, inQuenchTemperature, inHoldTime, inFinalTemperature, inCoolingTime):
        self.mQuenchTemperature = inQuenchTemperature
        self.mHoldTime = inHoldTime
        self.mFinalTemperature = inFinalTemperature
        self.mCoolingTime = inCoolingTime

    def Evaluate(self, inTime ):
        if ( inTime > (self.mCoolingTime + self.mHoldTime )):
            return self.mFinalTemperature
        elif ( inTime > self.mHoldTime ):
            return self.mQuenchTemperature+(self.mFinalTemperature-self.mQuenchTemperature)/self.mCoolingTime*(inTime-self.mHoldTime)
        else:
            return self.mQuenchTemperature
        
    def SafeTimeStep(self, inTime, inTemperatureStep):
        if ( inTime > (self.mCoolingTime + self.mHoldTime )):
            return self.mHoldTime * 0.001
        elif ( inTime > self.mHoldTime ):
            return inTemperatureStep * self.mCoolingTime / abs(  self.mFinalTemperature - self.mQuenchTemperature )
        else:
            return 0


########################################################################################
#     FixedTemperature class: inherits from Temperature class
########################################################################################

class FixedTemperature(Temperature):
    def __init__(self, inTemperature):
        self.mTemperature = inTemperature

    def Evaluate(self, inTime):
        return self.mTemperature
    
    def SafeTimeStep(self, inTime, inTemperatureStep = 0.001 ):
        return 0


########################################################################################
#     ZigzagTemperature class: inherits from Temperature class
########################################################################################

class ZigzagTemperature(Temperature):
    def __init__(self, inInitialTemperature, inLowTemperature, inHighTemperature, 
                 inFinalTemperature, inInitialCoolingTime, inReheatingTime, 
                 inFinalCoolingTime, inLowHoldTime, inHighHoldTime):
        
        self.mInitialTemperature = inInitialTemperature
        self.mLowTemperature = inLowTemperature  
        self.mHighTemperature = inHighTemperature  
        self.mFinalTemperature = inFinalTemperature  
        self.mInitialCoolingTime = inInitialCoolingTime
        self.mReheatingTime = inReheatingTime
        self.mFinalCoolingTime = inFinalCoolingTime
        self.mLowHoldTime = inLowHoldTime
        self.mHighHoldTime = inHighHoldTime

    def Evaluate(self, inTime):
        if ( inTime < 0 ):
            return self.mInitialTemperature
        elif (inTime < self.mInitialCoolingTime):
            return self.mInitialTemperature + ( self.mLowTemperature - self.mInitialTemperature ) / self.mInitialCoolingTime * inTime
        elif (inTime < (self.mInitialCoolingTime + self.mLowHoldTime)):
            return self.mLowTemperature
        elif (inTime < (self.mInitialCoolingTime + self.mLowHoldTime + self.mReheatingTime )):
            return self.mLowTemperature + ( self.mHighTemperature - self.mLowTemperature ) / self.mReheatingTime * ( inTime - self.mInitialCoolingTime - self.mLowHoldTime )
        elif ( inTime < (self.mInitialCoolingTime + self.mLowHoldTime + 
              self.mReheatingTime + self.mHighHoldTime) ):
            return self.mHighTemperature
        elif (inTime < (self.mInitialCoolingTime + self.mLowHoldTime +
              self.mReheatingTime + self.mHighHoldTime + self.mFinalCoolingTime )):
            return self.mHighTemperature + ( self.mFinalTemperature - self.mHighTemperature ) / self.mFinalCoolingTime * ( inTime - self.mInitialCoolingTime - self.mLowHoldTime - self.mReheatingTime - self.mHighHoldTime )
        else:
            return self.mFinalTemperature
    
    def SafeTimeStep(self, inTime, inTemperatureStep = 0.001 ):
        if ( inTime < 0 ):
            return 0
        elif ( inTime < self.mInitialCoolingTime ):
            return inTemperatureStep * self.mInitialCoolingTime / abs( self.mLowTemperature - self.mInitialTemperature )
        elif ( inTime < (self.mInitialCoolingTime + self.mLowHoldTime )):
            return self.mLowHoldTime * 0.001
        elif ( inTime < (self.mInitialCoolingTime + self.mLowHoldTime + self.mReheatingTime )):
            return inTemperatureStep * self.mReheatingTime / abs( self.mHighTemperature - self.mLowTemperature )
        elif ( inTime < (self.mInitialCoolingTime + self.mLowHoldTime + self.mReheatingTime + self.mHighHoldTime )):
            return self.mHighHoldTime * 0.001
        elif ( inTime < (self.mInitialCoolingTime + self.mLowHoldTime + self.mReheatingTime + self.mHighHoldTime + self.mFinalCoolingTime )):
            return inTemperatureStep * self.mFinalCoolingTime / abs( self.mHighTemperature - self.mFinalTemperature )
        else:
            return 0
    

########################################################################################
#     FromFileTemperature class: inherits from Temperature class
########################################################################################

class FromFileTemperature(Temperature):
    def __init__(self, inInitialTemperature, inFileName):
        self.mInitialTemperature = inInitialTemperature
        self.mFileName = inFileName
        self.initialize()

    def Initialize(self):
        try:
            fhTempPath = open(self.mFileName, 'r')
        except IOError:
            print('cannot open', self.mFileName)
        else:
            lines_temp = fhTempPath.readlines()
            NumPoints = len(lines_temp)
            TimePoints = np.zeros(NumPoints + 1)
            TemperaturePoints = np.zeros(NumPoints + 1)
            TimePoints[0] = 0
            TemperaturePoints[0] = self.mInitialTemperature
           
            for i in range(NumPoints):
                 line_temp_splitted = (lines_temp[i]).split()
                 TimePoints[i + 1] = float(line_temp_splitted[0])
                 TemperaturePoints[i + 1] = float(line_temp_splitted[1])
    
            NumPoints += 1                          
            fhTempPath.close()
            
        self.mNumPoints = NumPoints
        self.mTimePoints = TimePoints
        self.mTemperaturePoints = TemperaturePoints
        
    def Evaluate(self, inTime):
        if ( inTime >= self.mTimePoints[-1] ):
            return self.mTemperaturePoints[-1]
         
        if ( inTime <= 0 ):
            return self.mInitialTemperature
         
        for i in range(self.mNumPoints):
            if ( inTime == self.mTimePoints[i] ):
                return self.mTemperaturePoints[i]
            elif ( inTime < self.mTimePoints[i] ):
                return self.mTemperaturePoints[i-1] + ( self.mTemperaturePoints[i] - self.mTemperaturePoints[i-1] ) / (self.mTimePoints[i] - self.mTimePoints[i-1] ) * ( inTime - self.mTimePoints[i-1] )
            
    def SafeTimeStep(self, inTime, inTemperatureStep = 0.001):
        if ( inTime >= self.mTimePoints[-1] ):
            return 0

        if ( inTime < 0 ):
            return 0

        for i in range(self.mNumPoints):
            if ( inTime == self.mTimePoints[i] ) and ( i != self.mNumPoints):
                if ( ( self.mTemperaturePoints[i] != self.mTemperaturePoints[i - 1] ) and 
                    ( self.mTemperaturePoints[i + 1] != self.mTemperaturePoints[i] ) ):
                    return inTemperatureStep * 0.5 * (abs((self.mTimePoints[i] - self.mTimePoints[i - 1]) / (self.mTemperaturePoints[i] - self.mTemperaturePoints[i - 1])) + abs((self.mTimePoints[i + 1] - self.mTimePoints[i]) / (self.mTemperaturePoints[i + 1] - self.mTemperaturePoints[i])))
                elif ( self.mTemperaturePoints[i] != self.mTemperaturePoints[i - 1] ):
                    return inTemperatureStep * abs((self.mTimePoints[i] - self.mTimePoints[i - 1]) / (self.mTemperaturePoints[i] - self.mTemperaturePoints[i - 1]))
                elif ( self.mTemperaturePoints[i + 1] != self.mTemperaturePoints[i] ):
                    return inTemperatureStep * abs((self.mTimePoints[i + 1] - self.mTimePoints[i]) / (self.mTemperaturePoints[i + 1] - self.mTemperaturePoints[i]))
                else:
                    return 0
            elif ( inTime < self.mTimePoints[i] ):
                if ( self.mTemperaturePoints[i] != self.mTemperaturePoints[i - 1] ):
                    return inTemperatureStep * abs((self.mTimePoints[i] - self.mTimePoints[i - 1] ) / (self.mTemperaturePoints[i] - self.mTemperaturePoints[i - 1]))
            else:
                return ( self.mTimePoints[i] - self.mTimePoints[i - 1] ) * 0.001

        return 0
