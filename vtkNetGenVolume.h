//  Copyright (c) Kitware, Inc.  See license.md for details.

#ifndef vtkNetGenVolume_h
#define vtkNetGenVolume_h

#include "vtkUnstructuredGridAlgorithm.h"

/**\brief Generate a tetrahedral volume mesh from a surface using NetGen.
  *
  */
class vtkNetGenVolume : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkNetGenVolume* New();
  vtkTypeMacro(vtkNetGenVolume, vtkUnstructuredGridAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Specify the maximum global mesh size allowed
   */
  vtkSetClampMacro(MaxGlobalMeshSize, double, 0.0 , VTK_DOUBLE_MAX);
  vtkGetMacro(MaxGlobalMeshSize, double);
  //@}

  //@{
  /**
   * Specify the mesh density: 0...1 (0 => coarse; 1 => fine)
   */
  vtkSetClampMacro(MeshDensity, double, 0.0 , 1.0);
  vtkGetMacro(MeshDensity, double);
  //@}

  //@{
  /**
   * Generate second-order volume elements
   */
  vtkSetClampMacro(SecondOrder, int, 0 , 1);
  vtkGetMacro(SecondOrder, int);
  //@}

  //@{
  /**
   * Specify number of optimize steps to use for 3-D mesh optimization
   */
  vtkSetClampMacro(NumberOfOptimizeSteps, int, 0 , VTK_INT_MAX);
  vtkGetMacro(NumberOfOptimizeSteps, int);
  //@}

protected:
  vtkNetGenVolume();
  ~vtkNetGenVolume() override = default;

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int RequestData(vtkInformation* request,
    vtkInformationVector** inInfo,
    vtkInformationVector* outInfo) override;

  double MaxGlobalMeshSize;
  double MeshDensity;
  int SecondOrder;
  int NumberOfOptimizeSteps;

private:
  vtkNetGenVolume(const vtkNetGenVolume&) = delete;
  void operator=(const vtkNetGenVolume&) = delete;
};

#endif // vtkNetGenVolume_h
