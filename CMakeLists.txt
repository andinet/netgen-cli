cmake_minimum_required(VERSION 3.8)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project (Netgen)

set(VTK_LIBRARIES "")

find_package(VTK COMPONENTS
    vtkCommonColor
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkFiltersCore
    vtkFiltersExtraction
    vtkFiltersGeneral
    vtkFiltersGeometry
    vtkFiltersModeling
    vtkFiltersPoints
    vtkFiltersSources
    vtkFiltersVerdict
    vtkIOExodus
    vtkIOGeometry
    vtkIOImage
    vtkIOInfovis
    vtkIOLegacy
    vtkIOPLY
    vtkIOXML
    vtkImagingCore
    QUIET
)

set(KIT_LIBS ${VTK_LIBRARIES})

find_package(Netgen QUIET REQUIRED)

add_executable(NetGenVolumeMeshing NetGenVolumeMeshing.cxx vtkNetGenVolume.h vtkNetGenVolume.cxx)

target_link_libraries(NetGenVolumeMeshing ${KIT_LIBS} nglib)
  



