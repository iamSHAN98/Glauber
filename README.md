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

## Usage

### Compilation

Build the event generator binary `Glauber` using
the following commands

```shell
git clone --recursive ttps://github.com/iamSHAN98/Glauber.git
cd Glauber
mkdir build && cd build
cmake ..
make -jN
```

*Replace* `N` *with desired number of threads for
parallel compilation*

### Event Generation

A [configuration](config.yaml) in YAML format is given as
an input to the event generator binary

```shell
./Glauber /path/to/config.yaml
```

The key-value pairs in the configuration file are
quite self-explanatory (refer to the [report](Report.pdf)
for details). For example, the key `Monte Carlo`
takes a boolean value to choose between Optical
and Monte Carlo Glauber calculations. Events are
stored in a HDF5 file. Refer to `SetOutputFormat`
method in [Glauber.cc](src/Glauber.cc), [GlauberMC.cc](src/GlauberMC.cc)
and [Nucleus.cc](src/Nucleus.cc) for the quantities stored
event-by-event or use `h5dump` to list them from
the output file

```shell
h5dump -H /path/to/output
``` 

### Parallel Run

Large number of events can be generated using the
the [parallel script](parallel) provided

`./parallel --binary /path/to/Glauber --config /path/to/config.yaml --cores`*`N`*`--name`*`Directory`*

This initiates *`N`* parallel event generation
instances and saves the outputs in *`Directory`*
(will be created) as `Glauber_i.h5` with `i` from
1 to *`N`*. Use the optional `--input` to pass
additional arguments to the binary.

### [Analysis](tool/README.md)