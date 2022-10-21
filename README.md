### Making ntuples
The configuration file to use on MC: `TreeMaker/TM/test/treemaker_cfg_mc.py`.

To succesfully compile the ntuplizer code following steps are needed:
```
cmsrel CMSSW_12_4_0

cd CMSSW_12_4_0/src/

cmsenv

git cms-init

git clone git@github.com:pallabidas/TreeMaker.git

scram b -j 12
```
