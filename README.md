# Practical Cut Locus

This repository contains the implementation of the algorithm described in "Practical Computation of the Cut Locus on Discrete Surfaces ", C. Mancinelli, M. Livesu,  E. Puppo. The code consists of two apps:

1) one GUI (cutlocus_gui) supporting all the tools described in the paper,

2) one program (cutlocus_batch) that can be runned from command line and returns the cut locus as a .csv file.

## Compilation
You can compile the project using the standard CMake routine

`mkdir build`<br/>
`cd build`<br/>
`cmake ../`<br/>
`make`<br/>

For mac users, there is also the possibility of building the Xcode project by compiling a python script. Namely

`python3 scripts/build.py xcode`

This will generate the Xcode project "cutlocus" in the repository build/xcode.

## Dependencies
Mac users only need to have C++ 17 or later installed. Linux users may have to install some additional packages, as

`libxrandr-dev`<br/>
`libxinerama-dev`<br/>
`libxcursor-dev`<br/>
`libxi-dev.`<br/>

Similarly, Windows users may need to install some OS-related packages, 
but there are no dependencies that prevent the use of the code on such operating system. 

## Run
Once the projects are built, you can run both apps from within the project repository by issuing

`./bin/cutlocus_gui data/mesh.ply`

or

`./bin/cutlocus_batch data/mesh.ply vid VTP`

## Using the GUI

                                 ---Basics---
Rotate: left click + dragging.
Translate_ right click + dragging.                                 
                                 ---Source Selection---
There are two ways to select the source w.r.t which the cut locus must be computed. One way, is to
manually insert the index of the vertex in the dialog "Source", and then flag "Use Hardcoded Source".
Alternatively, a point on the mesh can be picked with shift + left_click. Even if the selected point is not a vertex of the mesh, the cut locus will be computed w.r.t the nearest vertex of the picked point.

                                ---Method for geodesic distances---
The method used to compute geodesic distances can be selected with the combobox "Method". The 
available choices are the one analyzed in the paper, i.e VTP, the Heat Method and a graph based 
geodesic solver.

                                ---Computing the cut locus---
Once the source has been chosen, press the button "Compute Cut Locus". The 
output will be the spanning tree described in the paper, pruned according to the values of the
enabled filter (see below). If the model has genus different from zero, by pressing "Close Loops" 
you will obtaining an approximation of the cut locus topologically consistent with the surface.
The button "Smoothing" allows you to smooth the result as described in the paper.

                                    ---Filtering---
Both the filtering techniques described in the paper can be enabled (and desabled) by flagging (and unflagging) the boxes "Grad Filter" and "Lap Filter". The policy for the filtering can be choose as well by acting on the radio buttons "Growing" and "Pruning". The thresholds can be set using the sliders "Teta" and "Lap". The extrema of the ranges in which such values vary can also be set acting on the corresponding dialogs ("Teta_Min","Teta_Max","Lap_Min" and "Lap_Max"). Similarly, you can remove branches shorter than the value reported by the slider "Depth" by flagging the checkbox "Cut Short Branches". The range of such slider can be set consistently with what seen for the gradient and the Laplacian. The "Thinning" check box is flagged by default, imposing the "zero-dimensionality" of the spanning tree (see the paper). Of course, it can be desactivated if needed.

                                 ---Smoothing parameters---
By default, everytime we perform smoothing we make 5 iterations following the gradient of the distance field, and 5 iterations following the discrete normal of the resulting polyline. The number of iterations of both procedure can be set in the section "Smoothing Parameters".

                                 ---Visualize---
Once the distance is computed, you can use the check box "Show Gradient" to visualize its gradient. There two sliders to set the size of the arrows and the lift along the normal of the mesh (useful when the surface has sharp or curvy features). 

## Using the bacth program

The batch program expect three arguments: the mesh, the index of a vertex of the mesh, and a method to compute the geodesic distances. For example, a proper input could be

`./bin/cutlocus_batch data/bunny_70k.ply 10 VTP`

Additionally, you can use some flags to enable all the tools described before. Namely
  

`options:`<br/>
  `--use_grad                    Use gradient filtering [false]`<br/>
  `--use_lap                     Use Laplacian filtering [false]`<br/>
  `--cut                         Cut short branches [false]`<br/>
  `--growing                     Set growing policy [false]`<br/>
  `--smoothed                    Smoothed CL [false]`<br/>
  `--use_grad_for_smoothing      Use gradient for smoothing [false]`<br/>

As you can see, all the options are disabled by default. Of course, if you enable a tool that require a value, you will have to write that value when asked from the program. For example, by running 

`./bin/cutlocus_batch data/bunny_70k.ply 10 VTP --use_grad`

you will be asked to enter the value of θ in the following way

Choose the threshold for the gradient filtering: <write the value and press RETURN>

The output will be one .csv "CL_On_Verts" containg the edges of the mesh forming of the cut locus (hence, pairs of indices of vertcies) and, if the flag --smoothed has been used, a .csv file "CL_Smoothed" with pairs of vec3f (the new positions of the vertices in "CL_On_Verts" after smoothing).




                         






