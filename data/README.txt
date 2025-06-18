alpha
beat
gamma
|B|: number of bodies
|M|: number of colors
|C|: number of configurations
|W|: number of options
|L|: number of lanes
q: per-lane capacity
c: color batch limit

-------------- d_bmj --------------- |B|*|M|*|C| rows
body: b	||	color: m	||	configuration: j	||	demand: d_bmj


-------------- sum d_bmj over on j ------------ |B|*|M| rows
body: b	||	color: m	||	demand: sum d_bmj over on j


-------------- r^w_j --------------- |C|*|W| rows
configuration: j		||	option: w		|| 	r^w_j


-------------- sigma_w ----------- |W| rows
option: w 		|| 		sigma_w
