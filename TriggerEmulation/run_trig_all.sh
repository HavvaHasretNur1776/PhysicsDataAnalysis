# lb-conda Semilep/hammer/2023-07-15 bash
#!/bin/bash

# Parameters you can change or pass as env vars or args if you want:
#years=("2016" "2017" "2018")
years=("2016")
magnets=("up")
BASE_PATH="/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/MCsamples/TrackerOnlyNTuples_Maryland/Angular_RDst_muonic_Run2/"

declare -A FILE_MAP=(
    ["signal"]="b0_dstmu/b0_dstmunu"
    #ADD REST OF SAMPLES
)

PYTHON_SCRIPT="run2-rdx-trg_emu.py"
#-> thi emulates all trig:hlt1+l0
#PYTHON_SCRIPT="run2-rdx-l0_global_tis.py"
#PYTHON_SCRIPT="run2-rdx-l0_hadron_tos.py"-> this havin problem in regression-looking for l0tosbranch- paradoxic
TREE="TupleB0/DecayTree"   # this is default
BMESON="b0"  # this is default
DEBUG_FLAG=""  # or "--debug"
LOAD='<triggers/l0/xgb4-2016.pickle>'
#LOAD='<triggers/l0/xgb4-2017.pickle>'
#LOAD='<triggers/l0/xgb4-2018.pickle>'
logs_dir="logs/"

mkdir -p "$logs_dir"

for year in "${years[@]}"; do
  for magnet in "${magnets[@]}"; do
    for key in "${!FILE_MAP[@]}"; do
      REL_PATH=${FILE_MAP[$key]}
      INPUT_FILE="${BASE_PATH}${year}/${REL_PATH}_${magnet}.root"
      OUTPUT_FILE="${BASE_PATH}${year}/${REL_PATH}_${magnet}_withAllTrigs.root"
      LOGFILE="${logs_dir}${year}_${key}_${magnet}.log"
      
      python3 $PYTHON_SCRIPT -t "$TREE" -y "$year" -B "$BMESON" -l "$LOAD" $DEBUG_FLAG "$INPUT_FILE" "$OUTPUT_FILE" > "$LOGFILE" 2>&1
      # For 17 and 18 specify change load files accordingly

      if [ $? -eq 0 ]; then
        echo "Successfully processed: $INPUT_FILE"
      else
        echo "Error processing: $INPUT_FILE (check $LOGFILE for details)"
      fi
    done
  done
done
