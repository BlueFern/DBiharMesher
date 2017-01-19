import vtk

ec_reader = vtk.vtkXMLPolyDataReader()
ec_reader.SetFileName("quadMeshFullECc4080.vtp")
ec_reader.Update()
ec_mesh = ec_reader.GetOutput()

atp_reader = vtk.vtkXMLPolyDataReader()
atp_reader.SetFileName("quadMeshFullATPc4080.vtp")
atp_reader.Update()
atp_map = atp_reader.GetOutput()
print atp_map.GetNumberOfCells()

atp_cell_data = atp_map.GetCellData().GetArray('ATP')

print atp_cell_data.GetTuple(0)

for i in range(320000,320000 + 80 * 80):
    previous_value = atp_cell_data.GetTuple(i - 3199)
    
    atp_cell_data.InsertNextTuple(previous_value)
    
print atp_cell_data.GetNumberOfTuples()
    
ec_mesh.GetCellData().AddArray(atp_cell_data)

writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName("quadMeshFullATPc4080_fixed.vtp")
writer.SetInputData(ec_mesh)
writer.Update()

