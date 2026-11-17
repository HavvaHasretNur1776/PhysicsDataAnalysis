import GaudiKernel.SystemOfUnits as Units
from Gaudi.Configuration import *
from Configurables import DaVinci, MCDecayTreeTuple

from Gaudi.Configuration import appendPostConfigAction
def doIt():
    """
    specific post-config action for (x)GEN-files
    """
    extension = "xgen"
    ext = extension.upper()

    from Configurables import DataOnDemandSvc
    dod  = DataOnDemandSvc ()
    from copy import deepcopy
    algs = deepcopy ( dod.AlgMap )
    bad  = set()
    for key in algs :
        if     0 <= key.find ( 'Rec'     )                  : bad.add ( key )
        elif   0 <= key.find ( 'Raw'     )                  : bad.add ( key )
        elif   0 <= key.find ( 'DAQ'     )                  : bad.add ( key )
        elif   0 <= key.find ( 'Trigger' )                  : bad.add ( key )
        elif   0 <= key.find ( 'Phys'    )                  : bad.add ( key )
        elif   0 <= key.find ( 'Prev/'   )                  : bad.add ( key )
        elif   0 <= key.find ( 'Next/'   )                  : bad.add ( key )
        elif   0 <= key.find ( '/MC/'    ) and 'GEN' == ext : bad.add ( key )

    for b in bad :
        del algs[b]

    dod.AlgMap = algs

    from Configurables import EventClockSvc, CondDB
    EventClockSvc ( EventTimeDecoder = "FakeEventTime" )
    CondDB  ( IgnoreHeartBeat = True )

appendPostConfigAction( doIt )

# 1. Basic DaVinci Settings
DaVinci().DataType   = "2016"
DaVinci().Simulation = True
DaVinci().EvtMax     = -1
DaVinci().TupleFile  = "MC-signal-2016-gen-level-only.root"

from Configurables import LoKi__Hybrid__MCTupleTool


mctuple = MCDecayTreeTuple("MCDecayTreeTuple")
# 1. Put the carets (^) DIRECTLY in the main decay string.
# This forces DaVinci to auto-generate the branches from the exact same parsed tree,
mctuple.Decay = "[B~0 -> ^(D*(2010)+ -> ^(D0 -> ^K- ^pi+) ^pi+) ^mu- ^nu_mu~]CC"

# 2. Add the tools globally.
# Because we used inline carets above, MCDecayTreeTuple will automatically 
# copy these tools down into EVERY auto-generated branch.
mctuple.ToolList = [
    "MCTupleToolHierarchy",
    "MCTupleToolKinematic",
    "LoKi::Hybrid::MCTupleTool/LoKi_MC"
]

# 3. Configure the LoKi variables once.
mctuple.addTool(LoKi__Hybrid__MCTupleTool, name="LoKi_MC")
mctuple.LoKi_MC.Variables = {
    "TRUE_PT"  : "MCPT",
    "TRUE_P"   : "MCP",
    "TRUE_ETA" : "MCETA"
}

# 6. Tell DaVinci to run our tuple
DaVinci().UserAlgorithms = [mctuple]

# 7. Load your .xgen file
from GaudiConf import IOHelper
IOHelper().inputFiles([
    './Gauss-11574021-1000ev-20260324.xgen' # Replace with your actual file name
], clear=True)