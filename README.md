# Glauber

Initial state event generator for heavy-ion collisions
based on [Optical and Monte Carlo Glauber](Report.pdf) model.

## Pre-requisites

### GSL

GSL is required for Gauss-Legendre quadrature used in
[Integration](include/Integration.h) module and sampling  distributions
used in [Random](src/Random.cc) module. The library and headers
can be installed using the Linux package manager

```shell
sudo apt install libgsl-dev
```

### YAML Parser

[Config](include/Config.h) module uses [yaml-cpp](https://github.com/jbeder/yaml-cpp) methods to
parse event generator configurations in YAML format.
yaml-cpp (added as Git submodule) is automatically
built during source compilation.

### DataStream

[DataStream](https://github.com/iamSHAN98/DataStream) is used to store generated events.
DataStream [pre-requisites](https://github.com/iamSHAN98/DataStream#installation) are required to be
installed first. DataStream (added as Git submodule)
is automatically built during source compilation.