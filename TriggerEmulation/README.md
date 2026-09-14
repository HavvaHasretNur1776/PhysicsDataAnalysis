# Emulation of L0Hadron, L0Global and HLT1 
This emulation strategy has been borrowed from Maryland github as following, refer to the [repo](https://github.com/umd-lhcb/TrackerOnlyEmu/tree/master) for further details:

```
git clone https://github.com/umd-lhcb/TrackerOnlyEmu.git
cd TrackerOnlyEmu
lb-conda Semilep/hammer/2023-07-15 bash # any python 3.8+ should work but NOT 3.12
make install # you can test but some reason L0 emulation not working but we dont need ...
cd scripts
cp ../run2-rdx-trg_emu.py #Use my version since there is incompatibility with ROOT6.28
cp ../run_trig_all.sh .
source run_trig_all.sh
```
## INFO
This script ```run2-rdx-trg_emu.py``` emulates ```d0_hlt1_trackmva_tos_emu```, ```d0_hlt1_twotrackmva_tos_emu``` through HLT1 step, ```d0_l0_hadron_tos_emu``` through XGB regression at L0hadron and ```b0_l0_global_tis_emu``` at L0Global step. See ```TrackerOnlyEmu/TrackerOnlyEmu/emulation/run2-rdx.py``` for the flow of these steps.

```
python3 run2-rdx-trg_emu.py -t "$TREE" -y "$year" -B "$BMESON" "$INPUT_FILE" "$OUTPUT_FILE"
```
### HLT1 Emulation process tree
- Apply Global Event Cut ([GEC](https://github.com/umd-lhcb/TrackerOnlyEmu/blob/master/TrackerOnlyEmu/triggers/hlt1/run2-Hlt1GEC.h)) to K, Pi
- [Emulate](https://github.com/umd-lhcb/TrackerOnlyEmu/blob/master/TrackerOnlyEmu/triggers/hlt1/run2-Hlt1TrackMVA.h) ```K_hlt1_trackmva``` and ```Pi_hlt1_trackmva```  
- Form ```d0_hlt1_trackmva``` -> ```K_hlt1``` OR ```Pi_hlt1```
- Form ```d0_hlt1_twotrackmva``` -> ```vector<track_spec[K,Pi]>``` (for the cuts refer [here](https://github.com/umd-lhcb/TrackerOnlyEmu/blob/master/TrackerOnlyEmu/triggers/hlt1/run2-Hlt1TwoTrackMVA.h))

### L0Global TIS emulation process
It is mimmicking the decision trhough efficiency file in ```TrackerOnlyEmu/TrackerOnlyEmu/triggers/l0/l0_tis_efficiency.root```. See ```triggers/l0/run2-L0Global.h``` for cut details.

### L0Hadron TOS emulation process
See ```TrackerOnlyEmu/TrackerOnlyEmu/triggers/l0/run2-L0Hadron.h``` for cut details.
XGB regression is used targetting this variable ```d0_l0_hadron_tos_emu```. There is already trained sample in ```TrackerOnlyEmu/TrackerOnlyEmu/triggers/l0/xgb4-{year}.pickle``` which is getting used in the main script to predict ```d0_l0_hadron_tos_emu```. There are other methods for regression like BDT but I didnt do training only used trained sample for predicting. See Maryland repo for further details. 
!!DO NOT FORGET SPECIFY SAMPLE NAME FOR 17 AND 18!!
