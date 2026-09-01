#1 1-dim each
### set model
drift = c("a*X+1","X+1")
# diffusion = matrix(c("0.3", "0", "0", "sigma"), nrow=2)
diffusion = matrix(c("b", "0", "0", "sigma"), nrow=2)
ymodel = setModel(
  drift=drift, 
  diffusion=diffusion, 
  solve.variable=c("X", "Y"),
  state.variable=c("X", "Y"),
  observed.variable="Y"
)

### set data
T <- 1
N <- 500
n <- N
h <- T/N
simulations_per_weight_update <- 1

true.par = list(
  a = -1.5,
  b = 0.7,
  sigma = 0.2
)
tmp.yuima <- simulate(ymodel, true.parameter=true.par, sampling=setSampling(n=N,Terminal=T))
ydata <- tmp.yuima@data

### set yuima
samp <- setSampling(delta = h, n=N)
variable_data_mapping <- list(
  "X" = 1,
  "Y" = 2
)


yuima <- setYuima(model = ymodel, data = ydata, sampling = samp, variable_data_mapping = variable_data_mapping)

n_particles = 1000
xinits <- matrix(numeric(n_particles), ncol = 1) # ncol = the number of unobserved variables, nrow = the number of particles

res <- simulate_multi_particles_with_weights(yuima, xinits=xinits, params = true.par, simulations_per_weight_update=simulations_per_weight_update)

filter_res <- kalmanBucyFilter(yuima,true.par,mean_init = 0,vcov_init = 0)

# plot weights against values
plot(res$values[,1,N+1],res$weights[,N+1])

# histogram of weights
bins <- hist(res$values[,1,N+1], plot = FALSE, breaks = 20)$breaks
grp <- cut(res$values[,1,N+1], breaks = bins, include.lowest = TRUE)
sum_w <- tapply(res$weights[,N+1], grp, sum)
labels <- names(sum_w)
barplot(
  sum_w,
  names.arg = labels,
  las = 2,
  xlab = "values",
  ylab = "sum of weights",
  main = "sum of weights",
  border = NA
)

#filter mean
filter_res@mean[N+1]

##approximated mean
approx_mean <- numeric(N+1)
for(i in 1:(N+1)){
  approx_mean[i] <- sum(res$values[,1,(i-1)*simulations_per_weight_update+1]*res$weights[,i])/sum(res$weights[,i])
}

plot((0:N)*h,approx_mean,type="l",col="red")
lines(filter_res@mean,col="blue")
lines(ydata@zoo.data$"Series 1")

# add legend, red: kalman, blue: particle, black: data, at outside of the plot
legend("topleft", legend=c("Particle", "Kalman", "Data"), col=c("red", "blue", "black"), lty=1, bty="n")