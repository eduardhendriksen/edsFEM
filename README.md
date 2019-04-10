## edsFEM

edsFEM is an open source, simple, 3d, frame analysis software package.

It is usable by setting up a structural model, including beams, supports and loads in a script.

Using these scripts it is possible to render model geometry and results of the analysis in a simple gui.

![3d gui](examples/example_industrial_support_3d.png)

## Features

- edsFEM supports static structural mechanical in three dimensions.
- A profile library containing all standard EU profiles is included in the profiles module.
- Renderer capable of rendering the model and results.
- View support conditions and member end releases in the gui.
- View applied loads and moments and the resulting moment, shear and normal force lines along the beams.
- Render and view displacements.
- Apply distributed loads across beams and surfaces easily.
- A simple optimizer to minimize structure weight.

## Examples

Examples of these scripts are found in the [examples](examples/) folder.

|3d Model industrial floor system| Loading |
|---|---|
|![3d floor system](examples/example_industrial_floor_3d.png)|![Loads](examples/example_industrial_floor_loads.png)|
|Moment lines about local y| Moment lines about local z|
|![My](examples/example_industrial_floor_My.png)|![Mz](examples/example_industrial_floor_Mz.png)|
|Normal stresses output| Displacement output |
|![My](examples/example_industrial_floor_sigma_n.png)|![Mz](examples/example_industrial_floor_u.png)|

## How it functions

edsFEM functions like a traditional FE-software. In its current iteration it limited to beam elements.

All user input, such as geometry, mechanical properties, supports and loads are assembled into a stiffness matrix and load vector for the entire system.

Using a solver supplied in the scipy package, the matrix equation containing the stiffness matrix, displacement vector and load vector is solved for the displacement vector.

Subsequently a postprocessor handles the process of making the results of solving the system viewable for the end user and determining internal forces and stresses in the system.

## How to install

To use edsFEM, you need to have the Python programming language installed.

To be able to use edsFEM to its full extents, you need the additional packages listed here;

|Functionality|packages|
|---|---|
|edsFEM|numpy|
||scipy|
|Renderer|panda3d|
||PIL|
