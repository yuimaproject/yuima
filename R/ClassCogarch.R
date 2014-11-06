# In order to deal the cogarch model in the \texttt{yuima} package 
# We define a new class called \texttt{yuima.cogarch} and its structure 
# is similar to those used for the carma model.
# The class \texttt{yuima.cogarch} extends the \texttt{yuima.model} and has
# an additional slot that contains informations about the model stored in an object of 
# class \texttt{cogarch.info}.
# The class \texttt{cogarch.info} is build  internally by the function \texttt{setCogarch} and it is a
# the first slot of an object of class \texttt{yuima.cogarch}.

# Class 'cogarch.info'
setClass("cogarch.info",
         representation(p="numeric",
                        q="numeric",
                        ar.par="character",
                        ma.par="character",
                        loc.par="character",
                        Cogarch.var="character",
                        V.var="character",
                        Latent.var="character",
                        XinExpr="logical",
                        measure="list",
                        measure.type="character")
)

# Class 'yuima.cogarch'

setClass("yuima.cogarch",
         representation(info="cogarch.info"),
         contains="yuima.model")

