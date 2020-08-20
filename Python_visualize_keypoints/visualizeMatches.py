#!/usr/bin/env python

import vtk
import os
import sys

# look here for various sources and their documentation
# http://www.vtk.org/doc/nightly/html/classvtkPolyDataAlgorithm.html
# geometry = "nested_sphere"

# Function definition is here
def make_vtk_model(filename_input, filename_output):
    count = 0
    scale = []
    color = []
    points = []
    r_prev = [0, 0, 0, 0]
    with open(filename_input) as f:
        for line in f:
            r = line.split('\t')
            print line + "\n"
            while len(r) > 3:
                try:
                    val = float(r[0])
                    print "value: " + str(val)
                    break
                except ValueError:
                    print "no good: " + r[0]
                    del r[0]
                    print "removed: " + r[0]

            if len(r) > 3:
                if r[1] <> r_prev[1] and r[2] <> r_prev[2] and r[0] <> r_prev[0]:
                    count = count + 1
                scale.append(float(r[3]))
                color.append(count)  # this could be disease status, etc
                points.append([-float(r[0]), float(r[1]), float(r[2])])
                r_prev = r

    print "Items read:", len(scale)

    # output objects
    outpd = vtk.vtkPolyData()
    outpoints = vtk.vtkPoints()
    outpointdata = outpd.GetPointData()

    out_array = vtk.vtkFloatArray()
    out_array.SetNumberOfComponents(1)
    out_array.SetName('ScaleFactor')
    print "Point data array created:", out_array.GetName(), out_array.GetNumberOfComponents()

    outpointdata.AddArray(out_array)
    # active array will cause scaling in glyph 3D
    outpointdata.SetActiveScalars('ScaleFactor')

    out_array2 = vtk.vtkFloatArray()
    out_array2.SetNumberOfComponents(1)
    out_array2.SetName('Color')
    print "Point data array created:", out_array2.GetName(), out_array2.GetNumberOfComponents()
    outpointdata.AddArray(out_array2)

    # fill in some data
    for s, p, c in zip(scale, points, color):
        if s < 15:
            out_array.InsertNextValue(s)
            out_array2.InsertNextValue(c)
            idx = outpoints.InsertNextPoint(p)
            print "Inserted value", out_array.GetValue(idx), out_array.GetTuple(idx), "at:", outpoints.GetPoint(idx)

    outpd.SetPoints(outpoints)

    # save it as visible spheres
    source1 = vtk.vtkSphereSource()
    source1.SetRadius(1)
    source1.SetPhiResolution(30)
    source1.SetThetaResolution(30)
    source1.Update()
    appender = vtk.vtkAppendPolyData()
    appender.AddInputData(source1.GetOutput())
    appender.Update()

    glyph = vtk.vtkGlyph3D()
    glyph.ScalingOn()
    glyph.SetInputData(outpd)

    glyph.SetSourceConnection(0, appender.GetOutputPort())
    glyph.Update()

    model_display = slicer.mrmlScene.CreateNodeByClass('vtkMRMLModelDisplayNode')
    slicer.mrmlScene.AddNode(model_display)
    print(model_display.GetID())
    model = slicer.mrmlScene.CreateNodeByClass('vtkMRMLModelNode')
    slicer.mrmlScene.AddNode(model)
    model.SetAndObserveDisplayNodeID(model_display.GetID())
    model.SetAndObservePolyData(source1.GetOutput())
    
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(glyph.GetOutput())

    writer.SetFileName(filename_output)
    writer.Write()

    return


if __name__ == "__main__":
    # Inputs
    filename_input = sys.argv[1]
    if not os.path.isfile(filename_input):
        print(str(filename_input) + " does not exists")
        exit(1)

    filename_output = filename_input.split(".txt")[0] + "_model.vtk"

    make_vtk_model(filename_input, filename_output)
 
