## We have splitted the simulate function into blocks to allow for future 
## methods to be added. S.M.I. & A.B.
## Interface to simulate() changed to match the S3 generic funciton in the 
## package stats
## added an environment to let R find the proper values defined in the main
## body of the function, which in turn calls different simulation methods
## All new simulation methods should look into the yuimaEnv for local variables
## when they need to "eval" R expressions

##:: function simulate
##:: solves SDE and returns result
subsampling <- function(x,y) return(x)

