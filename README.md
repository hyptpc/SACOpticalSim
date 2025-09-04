# How to build

```
mkdir build
cd build
cmake ..
make
```

# How to execute

after build

```
./SACOpticalSim <conf file> <output rootfile path> [macro]
```

for example

```
./SACOpticalSim ../conf/newSAC.conf test.root run.mac
```

when using BeamProfile, check the event number and set it in the mac file

```
[root] beam->GetEntries()
```