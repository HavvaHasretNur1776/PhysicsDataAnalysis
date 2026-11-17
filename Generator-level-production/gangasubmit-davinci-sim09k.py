import os

gauss_job_id = 0 
davinci_script_path = '/afs/cern.ch/work/h/havva/MC_Gen-level-prod/davinci_mc_options.py'

# 1. Build the official dataset
ds = LHCbDataset()
for sj in jobs(gauss_job_id).subjobs:
    if sj.status == 'completed':
        for f in sj.outputfiles.get(DiracFile):
            if f.namePattern.endswith('.xgen') and getattr(f, 'lfn', ''):
                ds.files.append(DiracFile(lfn=f.lfn))

if len(ds.files) == 0:
    print("Error: No LFNs found.")
    exit()

j = Job(name='2017_DaVinci_Grid_Final')

# 2. Use the modern Gaudi wrapper built for Ganga 8
app = prepareGaudiExec('DaVinci', 'v46r19')
app.platform = 'x86_64_v3-el9-gcc13-opt'
app.options = [davinci_script_path]

j.application = app
j.backend = Dirac()

# 3. This splitter ACTUALLY works in Ganga 8
j.inputdata = ds
j.splitter = SplitByFiles(filesPerJob=1)
j.outputfiles = [LocalFile('*.root')]

j.submit()