# crosshatching
tools for generating vector crosshatching drawings from raster images. Qt-based GUI, currently written to Visual Studio + Qt Tools extension. If there is interest I can change this to be CMake based. Contact me or file an issue if interested.

The codebase is also just me learning what can be done with the range-v3 library, which i've never touched before. Crosshatching is represented as lazy range views of polylines. All generation is done via views that are ultimately on top of iota/generate etc so that an entire "swatch" of crosshatching is an ephemeral object encoding a computation that will generate the swatch when iterated over.

![sample output](http://jwezorek.com/wp-content/uploads/2022/08/grace-drawing.png)
![sample output](http://jwezorek.com/wp-content/uploads/2022/08/castle-drawing.png)
