import os
import vtk
import argparse

OUT_DIR = '/hpc/scratch/bloodflo/Simulations/Flat/c4080f'
DATA_DIR = '/hpc/home/projects/bloodflo/simulations/CoupledCells3D/s.1.ATPv2.4K/solutionVTU'

ecFiles = [OUT_DIR + '/vtk/ec_mesh_parent.vtp', OUT_DIR + '/vtk/ec_mesh_left_daughter.vtp', OUT_DIR + '/vtk/ec_mesh_right_daughter.vtp']
smcFiles = [OUT_DIR + '/vtk/smc_mesh_parent.vtp', OUT_DIR + '/vtk/smc_mesh_left_daughter.vtp', OUT_DIR + '/vtk/smc_mesh_right_daughter.vtp']

ecAttribList = ['EC_Ca', 'EC_Ca_coupling', 'EC_IP3', 'EC_IP3_coupling', 'EC_SR', 'EC_Vm', 'EC_Vm_coupling']
smcAttribList = ['SMC_Ca', 'SMC_Ca_coupling', 'SMC_IP3', 'SMC_IP3_coupling', 'SMC_SR', 'SMC_Vm', 'SMC_Vm_coupling', 'SMC_w']

def mergeFiles(fileList):
    """Use vtkAppendPolyData filter to join data as read from the fileList.

    :param fileList: list of files to process.
"""
    appender = vtk.vtkAppendFilter()

    for file in fileList:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file)
        reader.Update()

        appender.AddInputData(reader.GetOutput())

    appender.Update()
   
    return appender.GetOutput()

if __name__ == '__main__':
    argParse = argparse.ArgumentParser(description='Copy attributes between geometries')
    argParse.add_argument('start', type=int, help='Start step')
    argParse.add_argument('end', type=int, help='End step')
    args = argParse.parse_args()

    ecMesh = mergeFiles(ecFiles)
    smcMesh = mergeFiles(smcFiles)
    
    if not os.path.exists('solution'):
        os.makedirs('solution')
    
    for step in range(args.start, args.end + 1):
    
        # Deal with EC data.
        ecFileName = '/ec_Data_t_' + str(step) + '.vtu'
        ecReader = vtk.vtkXMLUnstructuredGridReader()
        ecReader.SetFileName(DATA_DIR + ecFileName)
        ecReader.Update()
    
        ecData = ecReader.GetOutput()
    
        newEcMesh = vtk.vtkUnstructuredGrid()
        newEcMesh.DeepCopy(ecMesh)
    
        for attrib in ecAttribList:
            newEcMesh.GetCellData().AddArray(ecData.GetCellData().GetArray(attrib))
    
        print 'Wriing', ecFileName
        newEcWriter = vtk.vtkXMLUnstructuredGridWriter()
        newEcWriter.SetInputData(newEcMesh)
        newEcWriter.SetFileName(OUT_DIR + '/solution/' + ecFileName)
        newEcWriter.Update()
    
        # Deal with SMC data.
        smcFileName = '/smc_Data_t_' + str(step) + '.vtu'
        smcReader = vtk.vtkXMLUnstructuredGridReader()
        smcReader.SetFileName(DATA_DIR + smcFileName)
        smcReader.Update()
    
        smcData = smcReader.GetOutput()
    
        newSmcMesh = vtk.vtkUnstructuredGrid()
        newSmcMesh.DeepCopy(smcMesh)
    
        for attrib in smcAttribList:
            newSmcMesh.GetCellData().AddArray(smcData.GetCellData().GetArray(attrib))
    
        print 'Writing', smcFileName
        newSmcWriter = vtk.vtkXMLUnstructuredGridWriter()
        newSmcWriter.SetInputData(newSmcMesh)
        newSmcWriter.SetFileName(OUT_DIR + '/solution/' + smcFileName)
        newSmcWriter.Update()

