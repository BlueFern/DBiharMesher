import vtk

mapReader = vtk.vtkXMLPolyDataReader()
mapReader.SetFileName('ATPc4080.vtp')
mapReader.Update()

newMapArray = mapReader.GetOutput().GetCellData().GetArray('initialATP')

oldMapReader = vtk.vtkXMLPolyDataReader()
oldMapReader.SetFileName('quadMeshFullATPc4080.vtp')
oldMapReader.Update()

oldMap = oldMapReader.GetOutput()
oldMap.GetCellData().RemoveArray('initialATP')
oldMap.GetCellData().AddArray(newMapArray)

mapWriter = vtk.vtkXMLPolyDataWriter()
mapWriter.SetFileName('quadMeshFullATPc4080.vtp')
mapWriter.SetInput(oldMap)
mapWriter.Update()

