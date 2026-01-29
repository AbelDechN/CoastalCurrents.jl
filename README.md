# CoastalCurrents.jl

Contribution to the [BlueCloud 2026](https://blue-cloud.org/) project for the virtual lab [Coastal currents from observations](https://blue-cloud.d4science.org/group/coastalcurrentsfromobservations).
Here is a direct link to the [JupyterHub instance](https://jupyterhub.d4science.org/hub/oauth_login?context=%2Fd4science.research-infrastructures.eu%2FD4OS%2FCoastalCurrentsFromObservations).

You will need to have [Julia 1.9 or later](https://julialang.org/downloads/) installed. Start julia and issue the following command:


## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/AbelDechN/CoastalCurrents.jl")
```

```bash
pip install copernicusmarine
```


## Test download

Using a shell command:

```bash
# https://help.marine.copernicus.eu/en/articles/8632322-copernicus-marine-toolbox-troubleshoots
export COPERNICUSMARINE_CACHE_DIRECTORY=/tmp
copernicusmarine get --dataset-id cmems_obs-sl_eur_phy-ssh_my_al-l3-duacs_PT1S
```

Using julia:

```julia
ENV["COPERNICUSMARINE_CACHE_DIRECTORY"]="/tmp"
run(pipeline(`yes`,`copernicusmarine get --dataset-id cmems_obs-sl_eur_phy-ssh_my_al-l3-duacs_PT1S`))
```


