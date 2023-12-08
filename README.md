# Instructions on Set-up

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the git-history and may need to be downloaded independently.

2. Open a Julia console and do:

```
julia> using Pkg
julia> Pkg.activate("path/to/this/folder/on/your/computer")
julia> Pkg.instantiate()
```

This will install all necessary packages for you to be able to run the notebooks.

3. Open Jupyter notebook by doing the following:

```
julia> using IJulia
julia> notebook(; dir = ".")
```
