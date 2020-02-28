# Fibonacci encoding
## In short
Takes an input ROOT `TTree` with a branch: `std::vector<std::vector<short>> waveform` and compresses it. You can create such files by running Phil's module `WaveformToTreeDump_module` in `dunetpc/DAQSimAna`.

Modify the path of the `TFile` in the function `Compress` in `Compress.cxx`.

## To compile
```
mkdir build
cmake ..
make
```

## To run
```
./Compress
```
