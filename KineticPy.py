# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 12:38:05 2019

@author: Yihong Mauro
"""

from Simulation import *

if __name__ == "__main__":
    
    import time
    start_time = time.time()

    #theParameterFileName = "Se_test1.inp"
    theParameterFileName= input("Please enter the name of the input file:  ") 
    
    theSimulation = Simulation(theParameterFileName)    
    theSimulation.MasterEquationSolver()

    end_time = time.time()
    elapsed_time = time.time() - start_time
    print("The elapsed time for class-version program is %d seconds. \n"%(elapsed_time))
 
    