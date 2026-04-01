theApp.EvtMax = 1000

import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = ["EVNT.21991773._000007.pool.root.1","EVNT.21991773._000028.pool.root.1","EVNT.21991773._000031.pool.root.1","EVNT.21991773._000036.pool.root.1","EVNT.22035098._000024.pool.root.1","EVNT.22035098._000036.pool.root.1","EVNT.22035098._000058.pool.root.1","EVNT.22035321._000018.pool.root.1","EVNT.22035321._000029.pool.root.1","EVNT.22035321._000054.pool.root.1"]
#svcMgr.EventSelector.InputCollections = ["EVNT.22682785._000021.pool.root.1","EVNT.22682785._000023.pool.root.1","EVNT.22682785._000034.pool.root.1","EVNT.22682785._000035.pool.root.1","EVNT.22682785._000048.pool.root.1","EVNT.22682785._000053.pool.root.1","EVNT.22682785._000077.pool.root.1","EVNT.22682785._000085.pool.root.1","EVNT.22682785._000108.pool.root.1","EVNT.22682785._000117.pool.root.1"]
#svcMgr.EventSelector.InputCollections = ["EVNT.22682758._000001.pool.root.1","EVNT.22682758._000003.pool.root.1","EVNT.22682758._000004.pool.root.1","EVNT.22682758._000012.pool.root.1","EVNT.22682758._000037.pool.root.1","EVNT.22682758._000041.pool.root.1","EVNT.22682758._000072.pool.root.1","EVNT.22682758._000085.pool.root.1","EVNT.22682758._000109.pool.root.1","EVNT.22682758._000157.pool.root.1"]
#svcMgr.EventSelector.InputCollections = ["EVNT.14113765._000027.pool.root.1","EVNT.14113765._000319.pool.root.1","EVNT.14113765._000052.pool.root.1","EVNT.14113765._000388.pool.root.1","EVNT.14113765._000054.pool.root.1","EVNT.14113765._000399.pool.root.1","EVNT.14113765._000261.pool.root.1","EVNT.14113765._000402.pool.root.1","EVNT.14113765._000280.pool.root.1","EVNT.14113765._000478.pool.root.1"]
#svcMgr.EventSelector.InputCollections = ["EVNT.14113684._000048.pool.root.1","EVNT.14113684._000139.pool.root.1","EVNT.14113684._000178.pool.root.1","EVNT.14113684._000186.pool.root.1","EVNT.14113684._000212.pool.root.1","EVNT.14113684._000240.pool.root.1","EVNT.14113684._000247.pool.root.1","EVNT.14113684._000280.pool.root.1","EVNT.14113684._000294.pool.root.1","EVNT.14113684._000333.pool.root.1"]
#svcMgr.EventSelector.InputCollections = ["EVNT.22682838._000006.pool.root.1","EVNT.22682838._000014.pool.root.1","EVNT.22682838._000023.pool.root.1","EVNT.22682838._000028.pool.root.1","EVNT.22682838._000043.pool.root.1","EVNT.22682838._000100.pool.root.1","EVNT.22682838._000104.pool.root.1","EVNT.22682838._000108.pool.root.1","EVNT.22682838._000112.pool.root.1","EVNT.22682838._000166.pool.root.1"]
#svcMgr.EventSelector.InputCollections = ["EVNT.22682812._000161.pool.root.1","EVNT.22682812._000162.pool.root.1","EVNT.22682812._000163.pool.root.1","EVNT.22682812._000164.pool.root.1","EVNT.22682812._000165.pool.root.1","EVNT.22682812._000166.pool.root.1","EVNT.22682812._000167.pool.root.1","EVNT.22682812._000168.pool.root.1","EVNT.22682812._000169.pool.root.1","EVNT.22682812._000170.pool.root.1"]
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from xAODEventInfoCnv.xAODEventInfoCnvConf import xAODMaker__EventInfoCnvAlg
job += xAODMaker__EventInfoCnvAlg()

from Rivet_i.Rivet_iConf import Rivet_i
rivet = Rivet_i()
import os
rivet.AnalysisPath = os.environ['PWD']

rivet.Analyses += ['Wyy']
rivet.RunName = ''
rivet.HistoFile = 'WithCutsIncl_Wyy.yoda.gz'
#rivet.CrossSection = 86.064
rivet.IgnoreBeamCheck = True
rivet.SkipWeights=True
job += rivet
