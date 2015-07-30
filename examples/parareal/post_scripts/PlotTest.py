__author__ = 's.terzi'

from Plotter import *

nprocs = [
          # 2,
          # 4,
          8
          ]

plotData = [
            RunTypes.SDC_Fine,
            RunTypes.SDC_Coarse,
            # RunTypes.PARA_CLASSIC,
            RunTypes.PARA_HYBRID_FULL,
            RunTypes.PFASST,
            ]

plotTypes = [
             # PlotTypes.SpeedUp,
             # PlotTypes.Efficiency,
             PlotTypes.AllIterationErrors,
             # PlotTypes.AllIterationResiduals,
             # PlotTypes.LastStepError,
             # PlotTypes.LastIterationError,
             # PlotTypes.Iter
             ]

plot(plotTypes, plotData, nprocs)