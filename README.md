## Quick start
Create conda environment using mamba (optionally use conda command)
```
mamba create -n fitness-landscape-graph python=3.11
mamba activate fitness-landscape-graph
pip install polars networkx numpy logomaker matplotlib scipy tqdm  
```

## Gephi Setting to generate graphs
illustration graph Gephi preview setting:
- node:  Opacity: 90
- node labels:
	- font proportional size: yes
	- Apple Braille 4 Plain
	- no box
- edges:
	- thickness 1
	- rescale weight yes: min 3.0 max 10.0
	- edge opacity normally 80, set to 50 for high concentration of azt as there so many edges
- edge arrow size: 3
- edge labels
	- Apple Braille 20 Plain

landscape graph Gephi setting:
- node color: is_peak
- node size: 15 to 50, exponential spline
- edge color deterministic
- global virtual node fitness = 7.5
- concentration specific virtual node fitness = 6
- layout
  - ForceAtlas 2: scaling 150.0, prevent overlap
  - rotate such that the higher peak align to the right side
  - Network Spiltter 3D: Z-Maximum 30

z level is calculated like this: 
$$z_{level}= round(\frac{z_{V} \times z_{levels}}{z_{max}})$$
, where $z_{V}$ is the node's $[z]$ value; $z_{levels}$ is the user-defined maximum number of levels, and $z_{max}$ is highest value in $[z]$ column.
