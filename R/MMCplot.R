


MMCplot <- function(s, x, value=max(x), eq = 1, result){
  ptm=list()
  for (j in 1:ncol(s)) {
    sl = fastDummies::dummy_cols(s[, j])
    sl = sl[, -1]
    tt = c()
    for (k in 1:max(s)) {
      s_l = Hmisc::Lag(sl[, k])
      res = nnet::multinom(s[s_l == 1, eq] ~ x[s_l == 1], model = TRUE)
      v = stats::coefficients(res)

      if(max(s)>2){
        num=0
        for(i in 1:(max(s)-1)){
          ci = exp(sum(v[i,]*c(1,value)))
          num = num + ci
        }

        i=0
        rp = c()
        for(i in 1:max(s)){

          if(i == 1){
            px = 1/(1+num)
          }else{
            px = (exp(sum(v[i-1,]*c(1,value))))/(1+num)
          }

          rp = c(rp, px)
        }

      }else{
          num = exp(sum(v*c(1,value)))
        rp = c()
        for(i in 1:max(s)){

          if(i == 1){
            px = 1/(1+num)
          }else{
            px = (exp(sum(v*c(1,value))))/(1+num)
          }

          rp = c(rp, px)
        }
      }

      tt = rbind(tt, rp)
    }

    rownames(tt) <- NULL
    ptm[[j]] = tt

  }

  tpm = 0
  for(l in 1:ncol(s)){
    p = as.numeric(result[[eq + eq - 1]][l, 1]) * ptm[[l]]

    tpm = tpm + p
  }
  requireNamespace("markovchain")
  tpm_x = methods::new("markovchain", transitionMatrix = tpm)
  set.seed(123)
  markovchain::plot(tpm_x)
}
