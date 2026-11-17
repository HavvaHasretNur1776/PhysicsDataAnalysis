import os
# 200k events for year 16
nSubJobs = 20
eventsPerJob = 10000

#to learn which available platforms are suitable for specific gauss version
#lb-run --list Gauss/v49r21
j = Job(name='2016_Gauss_v49r21_11574021')#sim09k
j.backend = Dirac()
j.outputfiles = [DiracFile('*.xgen'), DiracFile('*.xml')]

app = Executable()
app.exe = '/bin/bash'

args_list = []
# Ensure this points to your actual decfile options!
input_files = [LocalFile('/afs/cern.ch/work/h/havva/Gen/DecFiles/options/11574021.py')]

for i in range(nSubJobs):
    # --- THIS IS THE UPDATED BLOCK ---
    opts_file = f'seed_opts_{i}.py'
    with open(opts_file, 'w') as f:
        f.write("from Configurables import LHCbApp\n")
        f.write(f"LHCbApp().EvtMax = {eventsPerJob}\n")
        f.write("from Gauss.Configuration import GenInit\n")
        f.write("GaussGen = GenInit('GaussGen')\n")
        f.write("GaussGen.RunNumber = 11574021\n")
        f.write(f"GaussGen.FirstEventNumber = {(i * eventsPerJob) + 1}\n")
    # ---------------------------------
    
    input_files.append(LocalFile(opts_file))

    # We are now officially time-traveling back to Scientific Linux 6
    cmd = (
        "lb-run -c x86_64-slc6-gcc49-opt Gauss/v49r21 gaudirun.py "
        "\$GAUSSOPTS/Gauss-Job.py "
        "\$GAUSSOPTS/Gauss-2016.py "
        "\$GAUSSOPTS/GenStandAlone.py "
        "11574021.py "
        "\$LBPYTHIA8ROOT/options/Pythia8.py "
        f"{opts_file}"
    )
    args_list.append(['-c', cmd])

app.args = args_list[0]
j.application = app

j.splitter = ArgSplitter(args=args_list)
j.inputfiles = input_files

j.submit()
