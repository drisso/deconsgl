# R package `deconsgl`

This package implements a deconvolution method based on sparse group lasso.

# Installation

The recommended way to install the `decongsl` package is the following.

Note that since this package depends on `lsgl`, which is currently only available via Github, `lsgl` needs to be installed first.

```{r}
library(devtools)
install_github("nielsrhansen/lsgl")
install_github("drisso/decongsl")
```
