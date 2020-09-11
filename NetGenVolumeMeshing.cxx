//  Copyright (c) Kitware, Inc.  See license.md for details.

#include "vtkNetGenVolume.h"

#include "vtkCleanPolyData.h"
#include "vtkDataSetReader.h"
#include "vtkNew.h"
#include "vtkSuperquadricSource.h"
#include "vtkTriangleFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"

int main(int argc, char* argv[])
{
  int status = 0;
  vtkNew<vtkSuperquadricSource> source;
  vtkNew<vtkTriangleFilter> filter;
  vtkNew<vtkCleanPolyData> merger;
  vtkNew<vtkNetGenVolume> mesher;
  vtkNew<vtkUnstructuredGridWriter> writer;
  source->SetThetaResolution(4);
  source->SetPhiResolution(4);
  source->ToroidalOn();
  filter->SetInputConnection(source->GetOutputPort());
  merger->SetInputConnection(filter->GetOutputPort());
  mesher->SetInputConnection(merger->GetOutputPort());
  writer->SetInputConnection(mesher->GetOutputPort());
  writer->SetFileName("output.vtk");
  writer->Write();

  auto* mesh = mesher->GetOutput();
  if (!mesh)
  {
    std::cerr << "No output.\n";
    status = 1;
    return status;
  }
  if (argc <= 1)
  {
    // Verify some output mesh properties
    if (mesh->GetNumberOfCells() != 143)
    {
      std::cerr << "Expected 143 cells, got " << mesh->GetNumberOfCells() << ".\n";
      status = 1;
    }
    if (mesh->GetNumberOfPoints() != 50)
    {
      std::cerr << "Expected 50 points, got " << mesh->GetNumberOfPoints() << ".\n";
      status = 1;
    }
  }
  return status;
}
