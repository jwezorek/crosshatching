# crosshatching
tools for generating vector crosshatching drawings from raster images. currently still  experimenting with the code i.e. you can generate vector crosshatched drawings with what is checked in right now but I don't have a cmd line interface etc.

The codebase is also just me learning what can be done with the range-v3 library, which i've never touched before. Crosshatching is represented as lazy range views of polylines. All generation is done via views that are ultimately on top of iota/generate etc so that an entire "swatch" of crosshatching is an ephemeral object encoding a computation that will generate the swatch when iterated over.

![sample output](http://jwezorek.com/wp-content/uploads/2022/03/grace_ch.png)
