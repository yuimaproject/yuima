
toLatex.yuima <- function(object,...){
	mod <- NULL
	if(class(object)=="yuima.model") mod <- object
	if(class(object)=="yuima") mod <- object@model
 	n.eq <- mod@equation.number
	dr <- paste("\\left(\\begin{array}{c}\n")

	for(i in 1:n.eq){
		dr <- paste(dr, substr(mod@drift[i],2,nchar(mod@drift[i])-1), "\\\\ \n")
	}
	dr <- paste(dr, "\\end{array}\\right)", sprintf("d%s",mod@time.variable))

	
	df <- paste("\\left(\\begin{array}{c}\n")
	
	for(i in 1:n.eq){
	  df <- paste(df, paste(mod@diffusion[[i]],collapse="+"))
		
	 
	df <- paste(df, "\\\\ \n")
        }
	df <- paste(df, "\\end{array}\\right)")

	n.n <- mod@noise.number
	wn <- paste("\\left(\\begin{array}{c}\n")
	wn <- paste(wn, paste(sprintf("dW%s",1:n.n), sep="", collapse="\\\\ "))
	wn <- paste(wn, "\n \\end{array}\\right)")

	st <- paste("\\left(\\begin{array}{c}\n")
	st <- paste(st, paste(sprintf("d%s",mod@solve.variable), sep="", collapse="\\\\ "))
	st <- paste(st, "\n \\end{array}\\right)")

	mysymb <- c("*", "alpha", "beta", "gamma", "delta","rho","theta","sqrt")
	myrepl<- c(" \\cdot ", "\\alpha ", "\\beta ", "\\gamma ", "\\delta ","\\rho ", "\\theta ", "\\sqrt ")
	ns <- length(mysymb)

	for(i in 1:ns){
 		dr <- gsub( mysymb[i], myrepl[i] , dr ,fixed="TRUE")
		df <- gsub( mysymb[i], myrepl[i] , df ,fixed="TRUE")
	}

#	cat("$$\n")
#	cat(st)
#	cat("\n = \n")
#	cat(dr)
#	cat("\n +\n")
#	cat(df)
#	cat("'\n")
#	cat(wn)
#	cat("\n $$\n")
#	cat("$$\n")
#	cat(sprintf("%s(0)=%f,\\quad",mod@solve.variable,mod@xinit))
#	cat("\n $$\n")

	body <- paste("%%% Copy and paste the following output in your LaTeX file")
	body <- body <- c(body,paste("$$"))
	body <- c(body, paste(st))
	body <- c(body, paste(" = "))
	body <- c(body, paste(dr))
	body <- c(body, paste(" +"))
	body <- c(body, paste(df))
	body <- c(body, paste("'"))
	body <- c(body, paste(wn))
	body <- c(body, paste("$$"))
	body <- c(body, paste("$$"))
	body <- c(body, paste(sprintf("%s(0)=%f,\\quad",mod@solve.variable,mod@xinit)))
	body <- c(body, paste("$$"))
	
	structure(body, class="Latex")
}

toLatex.yuima.model <- toLatex.yuima 
