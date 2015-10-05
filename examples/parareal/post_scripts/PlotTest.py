__author__ = 's.terzi'

from ProcessStarter import RunTypes
from Plotter import *

nprocs = [
          # 0,
          # 2,
          # 4,
          # 8,
          16,
          32,
          64
          ]

plotData = [
            RunTypes.SDC_Fine,
            RunTypes.SDC_Coarse,
            RunTypes.PARA_CLASSIC,
            RunTypes.PARA_HYBRID_PARTIAL,
            RunTypes.PARA_HYBRID_FULL,
            # RunTypes.PFASST,
            ]

plotTypes = [
             # PlotTypes.SpeedUp,
             PlotTypes.SpeedUpAlpha,
             # PlotTypes.Efficiency,
             # PlotTypes.AllIterationErrors,
             # PlotTypes.AllIterationResiduals,
             # PlotTypes.LastStepError,
             # PlotTypes.LastStepResidual,
             # PlotTypes.LastIterationError,
             # PlotTypes.Iter,
             # PlotTypes.ConvergenceData,
             # PlotTypes.Timings
             ]

dirs = [
        # 'runs_001',
        'runs_002',
        'runs_003',
        'runs_004'
        ]
for dir in dirs:
    plot(plotTypes, plotData, nprocs, root=dir)

# dirs = [
#         'ConvergenceDataSpatialCoarsening',
#         'ConvergenceDataTimeSpatialCoarsening'
#         ]
# for dir in dirs:
#     plot([PlotTypes.ConvergenceData], plotData, nprocs, root=dir)