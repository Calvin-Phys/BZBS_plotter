# Plot Calculated Band Structure on Measured ARPES Data

This tool can fetch calculated band structures from [Topological Materials Database](https://www.topologicalquantumchemistry.com/#/), and plot them on any Matlab figures. It can help compare the measured and calculated band structures and provide more insights into the ARPES data. 

Now the tool can only deal with Orthorhombic lattices. 


Basic Workflow:

0. Open `BZBS_plotter.mlapp`
1. Input `ICSD` of the selected material, and click "Get".
    - The `ICSD` of the material can be found on [Topological Materials Database](https://www.topologicalquantumchemistry.com/#/).
2. Wait until the information of the material appears in the text area.
    - Check if they match the desired material.
    - Can also preview the unit cell or the Brillouin zone.
3. Input the beginning and ending high symmetry point names in "Kpoint 1" and "Kpoint 2". Select "soc" or "no-soc" to indicate whether the calculation should include the spin-orbit coupling.
    - Click "Plot Band" to plot the selected band structure.
    - Click "Plot Band on Current Axes" to plot on the current active Matlab figure.
    - Click "Delete Lines" to remove all lines from the current active Matlab figure.
4. Plot Brillouin Zone (for Orthorhombic lattices): select the unit cell directions of the vertical and horizontal axis
    - "\\" means to ignore this axis

Components:

1. `BZBS_plotter.mlapp`: The main GUI
2. `bz_GRAPH.m`: The band structure profile form the Topological Materials Database
3. `cif2json.py`: parse the .cif file by [SeeK-Path](https://www.materialscloud.org/work/tools/seekpath)
