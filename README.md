# crosshatching
tools for generating vector crosshatching images from raster images. currently just experimenting with the code i.e. you can't actually generate vector crosshatched dra3wings with what is checked in.

The codebase is also just me learning what can be done with the range-v3 library, which i've never touched before. Crosshatching is represented as lazy range views of polylines. All generation is done via views that are ultimately on top of iota/generate etc so that an entire "swatch" of crosshatching is an ephemeral object encoding a computation that will generate the swatch when iterated over.
