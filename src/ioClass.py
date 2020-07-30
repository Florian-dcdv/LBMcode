from numpy import *; from numpy.linalg import *
import matplotlib.pyplot as plt; from matplotlib import cm
import numpy
import os
from pylab import *
import vtk
from vtk.util import numpy_support
import pickle

class ioLBM:
    def output(self,case,localVars,time,rho,u): 
    
        # Write VTK output
        if case.outputF=='VTK':
            filenameVTK = "./VTKResults/t"+str(time)+".vtk"
            print("\t File:",filenameVTK,"is writing...")
            velocity = vtk.vtkDoubleArray()
            velocity.SetNumberOfComponents(3)              
            density = vtk.vtkDoubleArray()
            density.SetNumberOfComponents(1)
            if case.FlowConfiguration=='MultiPhase':
                temperature = vtk.vtkDoubleArray()
                temperature.SetNumberOfComponents(1)
                pressure = vtk.vtkDoubleArray()
                pressure.SetNumberOfComponents(1)
                fsc = vtk.vtkDoubleArray()
                fsc.SetNumberOfComponents(3) 
            for j in range(0,case.ny):
                for i in range(0,case.nx):
                    ui, uj, uk = u[0,i,j], u[1,i,j], 0                    
                    velocity.InsertNextTuple3(ui, uj, uk)
                    density.InsertNextTuple([rho[i,j]])
                    if case.FlowConfiguration=='MultiPhase':
                        temperature.InsertNextTuple([localVars.Tim[i,j]])
                        pressure.InsertNextTuple([localVars.Peos[i,j]])
                        fsci, fscj, fsck = localVars.fsc[0,i,j], localVars.fsc[1,i,j], 0					
                        fsc.InsertNextTuple3(fsci, fscj, fsck)
            velocity.SetName('velocity')        
            density.SetName('rho')
            if case.FlowConfiguration=='MultiPhase':
                temperature.SetName('Temp')
                pressure.SetName('p')
                fsc.SetName('fsc')
            #File arrays and write into a VTK file 
            case.domain.GetPointData().SetVectors(velocity)                             
            case.domain.GetPointData().SetScalars(density)
            if case.FlowConfiguration=='MultiPhase':
                case.domain.GetPointData().AddArray(temperature)
                case.domain.GetPointData().AddArray(pressure)
                case.domain.GetPointData().AddArray(fsc)
            case.sg.SetFileName(filenameVTK)
            case.sg.SetInputData(case.domain)
            case.sg.Write()
        
        elif case.outputF=='screen':
            
            # Weird numbering of numpy array: use rot90
            #plt.clf(); plt.imshow(rho,cmap=cm.Reds)
            plt.clf(); plt.imshow(rot90(rho),cmap=cm.Reds)
            plt.colorbar()
            plt.title('Density')
            plt.show(block=False)
            if case.saveSnapshot==True:
                plt.savefig("./Snapshots/density_t"+str(time)+".png")
            plt.pause(2)
            plt.close()            
            
            plt.clf(); plt.imshow(rot90(sqrt(u[0,:,:]**2+u[1,:,:]**2)),cmap=cm.Reds)
            plt.colorbar()
            plt.title('Velocity')
            plt.show(block=False)
            if case.saveSnapshot==True:
                plt.savefig("./Snapshots/velocity_t"+str(time)+".png")
            plt.pause(2)
            plt.close()
            
            if case.FlowConfiguration=='MultiPhase':
                plt.clf(); plt.imshow(rot90(localVars.Tim),cmap=cm.Reds)
                plt.colorbar()
                plt.title('Temperature')
                plt.show(block=False)
                if case.saveSnapshot==True:
                    plt.savefig("./Snapshots/temperature_t"+str(time)+".png") 
                plt.pause(2)
                plt.close()            

                plt.clf(); plt.imshow(rot90(localVars.Peos),cmap=cm.Reds)
                plt.colorbar()
                plt.title('Peos')
                plt.show(block=False)
                if case.saveSnapshot==True:
                    plt.savefig("./Snapshots/Peos_t"+str(time)+".png")
                plt.pause(2)
                plt.close()             

        elif case.outputF=='line':
            # Open line output file
            filenameLine="./Lines/t"+str(time)+"_x="+str(case.lineXmin)+"-"+str(case.lineXmax) \
                                              +"y="+str(case.lineYmin)+"-"+str(case.lineYmax)+".dat"
            fline = open(filenameLine,"w+")
            print("\t File:",filenameLine,"is writing...")
            # Header
            if case.FlowConfiguration=='MultiPhase': 
                fline.write("# x \t y \t density \t velocity_x \t velocity_y \t temperature \t peos \n")
            else:
                fline.write("# x \t y \t density \t velocity_x \t velocity_y \n")
            for j in range(case.lineYmin,case.lineYmax+1):
                for i in range(case.lineXmin,case.lineXmax+1):
                    if case.FlowConfiguration=='MultiPhase': 
                        fline.write('{} \t {} \t {} \t {} \t {} \t {} \t {}\n'.format(i,j,rho[i,j],u[0,i,j],u[1,i,j],localVars.Tim[i,j],localVars.Peos[i,j]))
                    else:
                        fline.write('{} \t {} \t {} \t {} \t {} \n'.format(i,j,rho[i,j],u[0,i,j],u[1,i,j]))
                    #fline.write('{} \t {} \t {} \t {} \t {} \t {} \t {}\n'.format(j,i,rho[j,i],u[0,j,i],u[1,j,i],Tim[j,i],Peos[j,i]))
            fline.close() 
            
    def restart(self,globalVars,case,localVars,time):
        filenameRestart="t"+str(time)+".dat"
        print("\t File:","./Restarts/"+filenameRestart,"is writing...")
        output = open("./Restarts/"+filenameRestart, "wb") 
        if  case.solveTemp==True:
            timeFin = [time,localVars.fin,localVars.Tim]
        else:
            timeFin = [time,localVars.fin]    
        pickle.dump(timeFin,output)
        output.close()
        
                            
