__author__ = 's.terzi'

from Plotter import *
from ProcessStarter import RunTypes

nprocs = [
          # 0,
          # 2,
          # 4,
          # 8,
          16,
          # 32
          ]

plotData = [
            # RunTypes.SDC_Fine,
            # RunTypes.SDC_Coarse,
            RunTypes.PARA_CLASSIC,
            # RunTypes.PARA_HYBRID_FULL,
            RunTypes.PARA_HYBRID_PARTIAL,
            # RunTypes.PFASST,
            ]

plotTypes = [
             PlotTypes.SpeedUp,
             # PlotTypes.Efficiency,
             PlotTypes.AllIterationErrors,
             # PlotTypes.AllIterationResiduals,
             # PlotTypes.LastStepError,
             PlotTypes.LastIterationError,
             # PlotTypes.Iter
             ]

plot(plotTypes, plotData, nprocs)