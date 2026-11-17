## Generator level production

1. Create generator level .xgen file by running Gauss using custom decfile where I removed all generator-level cuts.


```
git lb-clone-pkg Gen/DecFiles -b Sim09
cd Gen/DecFiles
cp [your-dec-file].dec .
make
```
This will procude 11574021.py options file to be used at Gauss run. Ready Gauss env
```
lb-run Gauss/v49r21 bash
export DECFILESROOT=${PWD}
```

Then run Gauss (locally),

```
gaudirun.py $GAUSSOPTS/Gauss-Job.py $GAUSSOPTS/Gauss-2016.py $GAUSSOPTS/GenStandAlone.py \
            $DECFILESROOT/options/11574021.py $LBPYTHIA8ROOT/options/Pythia8.py myOptions.py
```

2. Run DaVinci over .xgen file to get the tuple 

```
lb-run -c best DaVinci/v46r19 gaudirun.py tuple-prod-signal-gen-level.py
```
