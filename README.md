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
./SACOpticalSim ../conf/newSAC.conf test.root vis.mac
```

newSAC.conf is for the new SAC (PMT 14ch) setup, and oldSAC.conf is for the old SAC (PMT 8ch).
