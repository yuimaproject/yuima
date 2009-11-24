# Class Definitions
# This source MUST be loaded first

# Class 'yuima.pars'

# parameter object included in 'yuima.model'
setClass("model.parameter",representation(all="character",
                                          common="character",
                                          diffusion="character",
                                          drift="character",
                                          jump="character",
                                          measure="character"
                                          )
         )

# Class 'yuima.model'
setClass("yuima.model",representation(drift="expression",
                                      diffusion="list",
                                      hurst="numeric",
                                      jump.coeff="expression",
                                      measure="list",
                                      measure.type="character",
                                      parameter="model.parameter",
                                      state.variable="character",
                                      jump.variable="character",
                                      time.variable="character",
                                      noise.number="numeric",
                                      equation.number="numeric",
                                      dimension="numeric",
                                      solve.variable="character",
                                      xinit="numeric"
                                      )
         )

# Class 'yuima.data'

# we want yuimaS4 to use any class of data as input
# the original data will be stored in OrigData
# we convert these objects internally to "zoo" object
# in the future, we may want to use more flexible
# classes

setClass("yuima.data", representation(original.data = "ANY",
                                      zoo.data = "ANY"
                                      )
         )


# Class 'yuima.sampling'

# sampling is now empty, but should give informations on the sampling
# type, rate, deltas, etc.

setClass("yuima.sampling", representation(Terminal = "numeric",
                                          division = "numeric",
										  Initial  = "numeric",
										  grid     = "numeric",
										  random   = "logical"
                                          )
         )

# Class 'yuima'

# this is the principal class of yuima project. It may contain up to
# three slots for now: the data, the model and the sampling

setClass("yuima.characteristic", representation(equation.number = "numeric",
                                                time.scale = "numeric"
                                                )
         )


setClass("yuima", representation(data = "yuima.data",
                                 model = "yuima.model",
                                 sampling = "yuima.sampling",
                                 characteristic = "yuima.characteristic"
                                 )
         )
