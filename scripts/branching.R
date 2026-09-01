# 各粒子は c(value, weight) という長さ2の数値ベクトルとする

branch <- function(particles) {
  n <- length(particles)
  if (n <= 1) {
    warning("number of particles should be greater than 1")
    return(particles)
  }
  us <- runif(n)
  new_particles <- list() # 新たな粒子はリストとして保持
  g <- n
  h <- n
  for (i in 1:(n - 1)) {
    g_int_part <- trunc(g)
    g_dec_part <- g - g_int_part
    nw <- particles[[i]][2] * n
    nw_int_part <- trunc(nw)
    nw_dec_part <- nw - nw_int_part
    u <- us[i]

    if (nw_dec_part <= g_dec_part) {
      if (u <= nw_dec_part / g_dec_part) {
        o <- nw_int_part + h - g_int_part
      } else {
        o <- nw_int_part
      }
    } else {
      if (u <= (1 - nw_dec_part) / (1 - g_dec_part)) {
        o <- nw_int_part + h - g_int_part
      } else {
        o <- nw_int_part + 1
      }
    }

    if (o > 0) {
      # 新たな粒子は元の粒子の value と 1/n の weight を持つ
      new_particle <- c(particles[[i]][1], 1 / n)
      # 同じ粒子を o 個複製して追加
      new_particles <- c(new_particles, replicate(o, new_particle,
        simplify = FALSE
      ))
    }
    g <- g - nw
    h <- h - o
  }
  if (h > 0) {
    new_particle <- c(particles[[n]][1], 1 / n)
    new_particles <- c(new_particles, replicate(h, new_particle,
      simplify = FALSE
    ))
  }
  new_particles
}

# テスト用関数
test <- function(n, debug = FALSE) {
  particles <- list()
  for (i in 1:n) {
    particles[[i]] <- c(i, runif(1))
  }
  # 重みの合計が 1 になるように正規化
  sum_weights <- sum(sapply(particles, function(x) x[2]))
  for (i in 1:n) {
    particles[[i]][2] <- particles[[i]][2] / sum_weights
  }
  print(paste("particles:", length(particles)))
  if (debug) {
    for (p in particles) {
      print(paste("value:", p[1], ", weight:", p[2], ", nw:", p[2] *
        length(particles)))
    }
  }
  new_particles <- branch(particles)
  print(paste("new_particles:", length(new_particles)))
  if (debug) {
    for (p in new_particles) {
      print(paste("value:", p[1], ", weight:", p[2], ", nw:", p[2] *
        length(new_particles)))
    }
  }
}

# 実行時間を計測する関数
measure_time <- function(n) {
  particles <- list()
  for (i in 1:n) {
    particles[[i]] <- c(i, runif(1))
  }
  sum_weights <- sum(sapply(particles, function(x) x[2]))
  for (i in 1:n) {
    particles[[i]][2] <- particles[[i]][2] / sum_weights
  }
  system.time(branch(particles))
}
