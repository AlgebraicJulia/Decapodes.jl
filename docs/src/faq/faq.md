# FAQ

## 1. How to incorporate scalar or vector field input data where you have a function of the embedded coordinates?

## 2. How to incorporate input data from a file with linear interpolation?

The Grigoriev Ice Cap model has a section where after the initial data is loaded from a TIF, the data is interpolated so that it may fit over a discrete mesh of our choosing. The link for that is [here](../grigoriev/grigoriev.md#loading-a-scientific-dataset).

## 3. How to set boundary conditions like fixed value, no flux, and no slip?

## 4. How to plot a derived quantity?

## 5. How to add artificial diffusion for 0- or 1-forms to improve stability?

## 6. How to use a Laplacian solver / multigrid?

To use multigrid methods in the Laplacian solver, you need to create a `PrimalGeometricMapSeries` that will take a coarse mesh and apply a subdivision method to it some number of times. After that, just use this result as you would a regular mesh for simulation.

```julia
  s = triangulated_grid(100,100,1,1,Point3D)

  series = PrimalGeometricMapSeries(s, binary_subdivision_map, 4);

  ...

  f_mg = sim_mg(series, generate);
```

## 7. How to do a bunch of workflows?
