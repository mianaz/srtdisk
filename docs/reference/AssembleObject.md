# Assemble an object from an h5Seurat file

Assemble an object from an h5Seurat file

## Usage

``` r
AssembleAssay(assay, file, slots = NULL, verbose = TRUE)

AssembleDimReduc(reduction, file, verbose = TRUE)

AssembleGraph(graph, file, verbose = TRUE)

AssembleImage(image, file, verbose = TRUE)

AssembleNeighbor(neighbor, file, verbose = TRUE)

AssembleSeuratCommand(cmd, file, verbose = TRUE)
```

## Arguments

- assay, reduction, graph, image, neighbor, cmd:

  Name of assay, reduction, graph, image, neighbor, or command to load

- file:

  A connected h5Seurat file to pull the data from

- slots:

  Optional vector of assay slots to load, defaults to all slots present
  in assay

- verbose:

  Show progress updates

## Value

`AssembleAssay`: An `Assay` object
