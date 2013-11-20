
toLatex.yuima <- function (object, ...) 
{
    mod <- NULL
    if (class(object) == "yuima.model") 
	mod <- object
    if (class(object) == "yuima") 
	mod <- object@model
    n.eq <- mod@equation.number
    dr <- paste("\\left(\\begin{array}{c}\n")
    for (i in 1:n.eq) {
        dr <- paste(dr, substr(mod@drift[i], 2, nchar(mod@drift[i]) - 
							   1), "\\\\ \n")
    }
    #
    dr <- paste(dr, "\\end{array}\\right)", sprintf("d%s", mod@time.variable))
    n.n <- mod@noise.number
    df <- paste(sprintf("\\left[\\begin{array}{%s}\n",paste(rep("c",n.n),sep="",collapse="")))
    for (i in 1:n.eq) {
        df <- paste(df, paste(mod@diffusion[[i]], collapse = "&"))
        df <- paste(df, "\\\\ \n")
    }
    df <- paste(df, "\\end{array}\\right]")
# We consider the Jump 6/11
    if (length(mod@jump.coeff)>=1){
      dj<-paste("\\left(\\begin{array}{c}\n")
      for (i in 1:n.eq) {
        #       dj <- paste(dj, substr(mod@jump.coeff[i], 2, nchar(mod@jump.coeff[i]) - 
        #                                1), "\\\\ \n")
        dj <- paste(dj, mod@jump.coeff[i], "\\\\ \n")
        
      }
      dj <- paste(dj, "\\end{array}\\right)", sprintf("d %s", mod@jump.variable))
    }
    
    wn <- paste("\\left(\\begin{array}{c}\n")
    wn <- paste(wn, paste(sprintf("dW%s", 1:n.n), sep = "", collapse = "\\\\ "))
    wn <- paste(wn, "\n \\end{array}\\right)")
    st <- paste("\\left(\\begin{array}{c}\n")
    st <- paste(st, paste(sprintf("d%s", mod@solve.variable), 
						  sep = "", collapse = "\\\\ "))
    st <- paste(st, "\n \\end{array}\\right)")
    mysymb <- c("*", "alpha", "beta", "gamma", "delta", "rho", 
				"theta","sigma","mu", "sqrt")
#     myrepl <- c(" \\cdot ", "\\alpha ", "\\beta ", "\\gamma ", 
# 				"\\delta ", "\\rho ", "\\theta ", "\\sqrt ")
    myrepl <- c(" \\cdot ", "\\alpha ", "\\beta ", "\\gamma ", 
                "\\delta ", "\\rho ", "\\theta ","\\sigma","\\mu", "\\sqrt ")
    ns <- length(mysymb)
    for (i in 1:ns) {
        dr <- gsub(mysymb[i], myrepl[i], dr, fixed = "TRUE")
        df <- gsub(mysymb[i], myrepl[i], df, fixed = "TRUE")
# for Jump         
        if (length(mod@jump.coeff)>=1){
          dj <- gsub(mysymb[i], myrepl[i], dj, fixed = "TRUE")
        }
    }
    body <- paste("%%% Copy and paste the following output in your LaTeX file")
    body <- body <- c(body, paste("$$"))
    body <- c(body, paste(st))
    body <- c(body, paste(" = "))
    body <- c(body, paste(dr))
    body <- c(body, paste(" +"))
    body <- c(body, paste(df))
    body <- c(body, paste("'"))
    body <- c(body, paste(wn))
    # For Jump 6/11
    if (length(mod@jump.coeff)>=1){
      body <- c(body, paste(" +"))
      body <- c(body, paste(dj))
    }
    
    body <- c(body, paste("$$"))

    body <- c(body, paste("$$"))
# For Initial Conditions     
    numb.solve.var <- length(mod@solve.variable)
    bodyaus <-c( paste("\\left(\\begin{array}{c}\n"))
    for (i in 1:numb.solve.var) {
      bodyaus <- paste(bodyaus, paste(paste(mod@solve.variable[i],"(0)",sep=""),substr(mod@xinit[i], 2, nchar(mod@xinit[i]) - 
                                                                                         1),sep="="), "\\\\ \n")
    }
    bodyaus <- paste(bodyaus, "\\end{array}\\right)")
    for (i in 1:ns) {
      bodyaus <- gsub(mysymb[i], myrepl[i], bodyaus, fixed = "TRUE")
    }
    
    body<-c(body,paste(bodyaus))
    
#     body <- c(body, paste(sprintf("%s(0)=%f,\\quad", mod@solve.variable, 
# 								  mod@xinit)))
    body <- c(body, paste("$$"))
    structure(body, class = "Latex")
}



toLatex.yuima.model <- toLatex.yuima 
