## edsFEM

edsFEM is an open source, simple, 3d, frame analysis software package.

It is usable by setting up a structural model, including beams, supports and loads in a script.

Using these scripts it is possible to render model geometry and results of the analysis in a simple gui.

![3d gui](examples/example_industrial_support_3d.png)

Examples of these scripts are found in the [examples](examples/) folder.

|3d Model industrial floor system| Loading |
|---|---|
|![3d floor system](examples/example_industrial_floor_3d.png)|![Loads](examples/example_industrial_floor_loads.png)|
|Moment lines about local y| Moment lines about local z|
|![My](examples/example_industrial_floor_My.png)|![Mz](examples/example_industrial_floor_Mz.png)|
|Normal stresses output| Displacement output |
|![My](examples/example_industrial_floor_sigma_n.png)|![Mz](examples/example_industrial_floor_u.png)|

To use edsFEM, you need to have the Python programming language installed.

To be able to use edsFEM to its full extents, you need the additional packages listed here;

|Functionality|packages|
|---|---|
|edsFEM|numpy|
||scipy|
|Renderer|panda3d|
||PIL|
