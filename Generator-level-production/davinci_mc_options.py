import GaudiKernel.SystemOfUnits as Units
from Gaudi.Configuration import *
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

from Configurables import DaVinci, MCDecayTreeTuple
from DecayTreeTuple.Configuration import *

# =============================================================================
# 1. Configure the MC Truth Tuple
# =============================================================================
mctuple = MCDecayTreeTuple("MCDecayTreeTuple")

# The fully optimized decay descriptor:
# - Outer CC: Captures both B0 and B~0 decays.
# - Inner [ ]cc: Protects against B0 -> B~0 mixing drops.
# - Loose arrows (==>): Prevents EvtGen/Photos radiative photons from breaking the selection.
# - Carets (^): Automatically generates branches for the kinematics of these specific particles.
#mctuple.Decay = "[ [B~0]cc ==> ^(D*(2010)+ ==> ^(D0 ==> ^K- ^pi+) ^pi+) ^mu- ^nu_mu~ ]CC"
#mctuple.Decay = "[B~0 -> ^(D*(2010)+ -> ^(D0 -> ^K- ^pi+) ^pi+) ^mu- ^nu_mu~]CC"
# Mixing safe version: gaining 165 events in 1k
mctuple.Decay = "[ [B~0]cc -> ^(D*(2010)+ -> ^(D0 -> ^K- ^pi+) ^pi+) ^mu- ^nu_mu~]CC"

# Standard suite of tools for generator-level kinematic analysis
mctuple.ToolList = [
    "MCTupleToolKinematic",  # Gets P, PT, E, M, etc.
    "MCTupleToolHierarchy",  # Gets Mother/Daughter relationship info
    "MCTupleToolEventInfo"   # Gets the Run Number and Event Number
]

# =============================================================================
# 2. Configure DaVinci
# =============================================================================
DaVinci().Simulation = True
DaVinci().Lumi = False       # Lumi is meaningless for pure gen-level
DaVinci().DataType = "2016"  # Must match your Gauss generation year
DaVinci().UserAlgorithms = [mctuple]

# Process every single event in whatever .xgen file Ganga feeds it
DaVinci().EvtMax = -1  

# Name the output file (Ganga will automatically grab this and bring it back)
#DaVinci().TupleFile = "11574021_tuple.root" # get from Gauss job number --- IGNORE ---