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
git clone --recursive https://github.com/iamSHAN98/Glauber.git
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

```shell
./parallel --binary /path/to/Glauber --config /path/to/config.yaml --cores N --name Directory
```

This initiates N parallel event generation
instances and saves the outputs in *Directory*
(will be created) as `Glauber_i.h5` with `i` from
1 to N. Use optional argument `--input` to pass
additional parameter(s) to the binary.

### [Analysis](tool/README.md)

## Extension

### Add Nucleus

<details>
<summary> Expand </summary>

`Nucleus` key in configuration file takes `string`
as inputs that are mapped to C `enum KeyNucleus`
using STL Map `MapNucleus` (refer to [KeyMap.h](include/KeyMap.h)).
These `enum` variables are used as indices to pick
the corresponding nuclear parameters from an array of
`struct NucleusData` (refer to [NucleusData.h](include/NucleusData.h)).
Nucleon density profiles are defined as methods of
[ChargeProfile](include/Profile.h). Given a `NucleusData` instance,
corresponding profile is chosen based on C `enum 
KeyChargeProfile` from a `switch-case` (refer to
[Profile.cc](src/Profile.cc)). There are two possible scenarios
when adding a new nucleus,

* For an existing nucleon density profile, just add
  an entry for `KeyNucleus`, `MapNucleus` in KeyMap.h
  and a `NucleusData` instance in NucleusData.h.

* In case of a new profile, check if defintiion of 
  `NucleusData` needs to be changed or not. Then,

  1. Define the corresponding method in `ChargeProfile`.
  2. Add an entry for `KeyChargeProfile` in KeyMap.h
  	 and a `case` statement in `GetProfile`.

  Now, instructions for the first scenario are to
  be followed only.

</details>

### Use As Module

<details>
<summary> Expand </summary>

There are two possible use cases in this regard,

* Many initial state models use Glauber calculations
  as a stepping stone - utilizes the thickness function
  or sampled nucleon positions. Use either `Glauber`
  or `GlauberMC` as base class in this case. Apart
  from defining model-specific new methods in the
  derived class, base class methods (e.g. `GenerateEvent`)
  may need to be overridden. To be noted, `GlauberMC`
  itself inherits from `Glauber`.

* Glauber model can be used to generate initial
  condition for hydrodynamic simulation to model
  space-time evolution of Quark Gluon Plasma (QGP)
  produced in heavy-ion collisions. Use method
  `GetTransverseProfile` for each node of a numerical
  grid in this case. Grid specs for the transverse
  plane can be passed as command line arguments to
  the binary. Accordingly modify [main](main.cc),

  ```cpp
  while(n < Nev){
    if(Obj->GenerateEvent()){
      .
      .
      .
      EventGen::Position R;
        
      for(int i = 0; i < SizeX; i++){
        R.x = MinX + i*StepX;

        for(int j = 0; j < SizeY; j++){
          R.y = MinY + j*StepY;

          Profile[i][j] = Obj->GetTransverseProfile(R);
          .
          .
          .
  ```

</details>

### [Add Plotting Script](tool/README.md#Add-New-Script)
