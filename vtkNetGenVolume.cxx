//  Copyright (c) Kitware, Inc.  See license.md for details.

#include "vtkNetGenVolume.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVector.h"
#include "vtkVectorOperators.h"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4067) /* unexpected tokens following preprocessor directive */
#endif
namespace nglib
{
#include "nglib.h"
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <array>
#include <map>

vtkStandardNewMacro(vtkNetGenVolume);

vtkNetGenVolume::vtkNetGenVolume()
{
  this->MaxGlobalMeshSize = 1.0e+6;
  this->MeshDensity = 0.5;
  this->SecondOrder = 0;
  this->NumberOfOptimizeSteps = 3;
}

void vtkNetGenVolume::PrintSelf(std::ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "MaxGlobalMeshSize: " << this->MaxGlobalMeshSize << "\n";
  os << indent << "MeshDensity: " << this->MeshDensity << "\n";
  os << indent << "SecondOrder: " << (this->SecondOrder ? "true" : "false") << "\n";
  os << indent << "NumberOfOptimizeSteps: " << this->NumberOfOptimizeSteps << "\n";
}

int vtkNetGenVolume::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
  {
    // We take polydata in (and don't override FillOutputPortInformation
    // so vtkUnstructuredGrids are output).
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  }
  else
  {
    return this->Superclass::FillInputPortInformation(port, info);
  }

  return 1;
}

int vtkNetGenVolume::RequestData(vtkInformation* /*request*/,
    vtkInformationVector** inInfo,
    vtkInformationVector* outInfo)
{
  int status = 0;
  auto* constraints = vtkPolyData::GetData(inInfo[0], 0);
  if (!constraints)
  {
    vtkErrorMacro("Nothing to mesh.");
    return status;
  }
  auto* points = constraints->GetPoints();

  auto* output = vtkUnstructuredGrid::GetData(outInfo, 0);
  if (!output)
  {
    vtkErrorMacro("No place to store result.");
    return status;
  }

  // Initialize netgen and create a mesh.
  nglib::Ng_Init();
  try
  {
    auto* mesh = nglib::Ng_NewMesh();

    vtkIdType numberOfPoints;
    vtkIdType const* connectivity;
    vtkVector3d coords;
    std::map<vtkIdType, int> pointMap;
    int nextPoint = 1;

    // Add vertex constraints
    vtkCellArray* verts = constraints->GetVerts();
    for (verts->InitTraversal(); verts->GetNextCell(numberOfPoints, connectivity); )
    {
      // vtkPolyData's vertices can be poly-vertices
      for (vtkIdType ii = 0; ii < numberOfPoints; ++ii)
      {
        if (pointMap.find(connectivity[ii]) == pointMap.end())
        {
          points->GetPoint(connectivity[ii], coords.GetData());
          nglib::Ng_AddPoint(mesh, coords.GetData());
          pointMap[connectivity[ii]] = nextPoint++;
        }
      }
    }

    // Add curve constraints (which we turn into point constraints)
    vtkCellArray* edges = constraints->GetLines();
    for (edges->InitTraversal(); edges->GetNextCell(numberOfPoints, connectivity); )
    {
      // vtkPolyData's lines can be poly-lines
      if (connectivity[0] == connectivity[numberOfPoints - 1])
      {
        // The curve is periodic... don't insert the same vertex twice
        --numberOfPoints;
      }
      for (vtkIdType ii = 0; ii < numberOfPoints; ++ii)
      {
        if (pointMap.find(connectivity[ii]) == pointMap.end())
        {
          points->GetPoint(connectivity[ii], coords.GetData());
          nglib::Ng_AddPoint(mesh, coords.GetData());
          pointMap[connectivity[ii]] = nextPoint++;
        }
      }
    }

    // Add surface constraints (ignore all but triangles and quads)
    vtkCellArray* polys = constraints->GetPolys();
    for (polys->InitTraversal(); polys->GetNextCell(numberOfPoints, connectivity); )
    {
      if (numberOfPoints < 3 || numberOfPoints > 4)
      {
        static bool once = false;
        if (!once)
        {
          once = true;
          vtkWarningMacro("Skipping all polygons except triangles and quads.");
        }
        continue;
      }
      // Ensure points have been added to mesh and translate connectivity into netgen point IDs.
      std::array<int, 4> ngconn;
      for (vtkIdType ii = 0; ii < numberOfPoints; ++ii)
      {
        if (pointMap.find(connectivity[ii]) == pointMap.end())
        {
          points->GetPoint(connectivity[ii], coords.GetData());
          nglib::Ng_AddPoint(mesh, coords.GetData());
          pointMap[connectivity[ii]] = nextPoint++;
        }
        ngconn[ii] = pointMap[connectivity[ii]];
      }
      // Insert a surface element
      nglib::Ng_AddSurfaceElement(mesh, numberOfPoints == 3 ? nglib::NG_TRIG : nglib::NG_QUAD, ngconn.data());
    }

    // Mesh the volume
    nglib::Ng_Meshing_Parameters mp;
    mp.maxh = this->MaxGlobalMeshSize;
    mp.fineness = this->MeshDensity;
    mp.second_order = this->SecondOrder;
    mp.optsteps_3d = this->NumberOfOptimizeSteps;
    nglib::Ng_Result result = nglib::Ng_GenerateVolumeMesh(mesh, &mp);

    if (result == nglib::NG_OK)
    {
      status = 1;
      // Add points to output
      points = output->GetPoints();
      if (!points)
      {
        points = vtkPoints::New();
        output->SetPoints(points);
      }
      numberOfPoints = nglib::Ng_GetNP(mesh);
      points->SetNumberOfPoints(numberOfPoints);
      for (int pp = 0; pp < numberOfPoints; ++pp)
      {
        nglib::Ng_GetPoint(mesh, pp + 1, coords.GetData());
        points->SetPoint(pp, coords.GetData());
      }

      // Add volume elements to output.
      std::array<int, 10> ngconn; // Current max number of points per element is 10
      std::array<vtkIdType, 10> vtkconn;
      int cellType;
      int nconn;
      int numberOfElements = nglib::Ng_GetNE(mesh);
      output->Allocate(numberOfElements);
      for (int ee = 0; ee < numberOfElements; ++ee)
      {
        auto elementType = nglib::Ng_GetVolumeElement(mesh, ee + 1, ngconn.data());
        switch (elementType)
        {
        case nglib::NG_TET:
          nconn = 4;
          cellType = VTK_TETRA;
          break;
        case nglib::NG_PYRAMID:
          nconn = 5;
          cellType = VTK_PYRAMID;
          break;
        case nglib::NG_PRISM:
          nconn = 6;
          cellType = VTK_WEDGE;
          break;
        case nglib::NG_TET10:
          nconn = 10;
          cellType = VTK_QUADRATIC_TETRA;
        default:
          vtkErrorMacro("Unknown output mesh element type (" << elementType << ").");
          break;
        }
        for (int ii = 0; ii < nconn; ++ii)
        {
          vtkconn[ii] = static_cast<vtkIdType>(ngconn[ii] - 1);
        }
        output->InsertNextCell(cellType, nconn, vtkconn.data());
      }
    }
    else
    {
      vtkErrorMacro("Could not generate volume mesh.");
    }
  }
  catch (...)
  {
    status = 0;
    vtkErrorMacro("Meshing failed.");
  }

  nglib::Ng_Exit();
  return status;
}
