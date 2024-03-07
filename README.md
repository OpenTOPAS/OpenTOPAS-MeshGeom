# Tetrahedral Mesh Geometry Capabilities for TOPAS

<p align="center">
  <a href="https://github.com/mghro/TOPAS-MeshGeom"><img alt="Adult Female ICRP 145 mesh head" src="https://github.com/mghro/TOPAS-MeshGeom/assets/34250777/18f99b31-6d5e-4f20-bf08-f3111a0a99f6" width="40%"></a>
  <p align="center"></p>
</p>

Largely based on the [Geant4 example](https://github.com/Geant4/geant4/tree/dd1f179cda58f54140945ad67846ff417903a862/examples/advanced/ICRP145_HumanPhantoms) using the ICRP 145 phantoms 



Authors:
 - [Isaac Meyer](imeyer@mgh.harvard.edu)
 - [Alejandro Bertolet](abertoletreina@mgh.harvard.edu)
 - [Dohyeon Yoo](dhyoo@yuhs.ac)

# Loading mesh geometries
Currently only the .ele/.node/.material format of tetrahedral meshes is supported where the file formats of ICRP 145 are followed.
## MRCP Phantoms
Phantoms need to be defined as geometry components in the parameter file (*Type="TsMRCP"*), and the directory where the files .node, .ele and .material provided by IRCP are located needs to be specified (see example "15F_SOBP_Liver.txt"). The mesh phantom geometries can be downloaded from the [ICRP website](https://www.icrp.org/publication.asp?id=ICRP%20Publication%20145). 
## Omed Phantom
Included in the `examples` folder is a coarse mesh geometry phantom named Omed (demo backwards). Visualization of the ICRP phantoms is extremely slow in the TOPAS GUI, and not recommended for checking geometry setups. Omed is oriented in the same direction as the ICRP 145 phantoms with a floating "R" for a right arm and "L" for left arm.
## Other Tetrahedral Meshes
Can be defined using the specification of file paths for each file e.g.

```
s:Ge/TetGeom/PhantomDirectory = "/mesh/file/directory/"
s:Ge/TetGeom/NodeFile         = "my_mesh.node"
s:Ge/TetGeom/MaterialFile     = "my_mesh.material"
s:Ge/TetGeom/EleFile          = "my_mesh.ele"
```

# Scoring
A specific class of scorer is needed for MRCP phantoms by specifying (*Quantity="TsMRCPScorer"*). By default, dose to water is computed, but other material can be used to calculate the dose by using the parameter *Material*.
For each organ involved it is necessary to add a new scorer, specifying the list parameter *Organ=1 "Name"*, where *Name* should be one (or more) of the materials listed in the .material file of the MRCP phantom, e.g., "Liver" or "Stomach_contents". This will restrict the dose considered only to the organ(s) of interest.
