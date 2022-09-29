


MMCXplot <- function(s, x, mov = TRUE, order=5, wdth=7, hght = 7, unts='in'){

out=c()
m1 = max(s)
t = nrow(s)
time = c()
px = c()

plots <- function(out){

  all = 1:length(out)
  return = paste('gridExtra::arrangeGrob( ', paste(paste('eval(as.name(out[',all,']))',sep=''), collapse = ', '), ')', sep='')

  return(return)
}

for (i in 1:ncol(s)){
 for (j in 1:ncol(s)){
    #create dummy variables of St-1
    sl = fastDummies::dummy_cols(s[, j])
    sl = sl[, -1]

    for (k in 1:m1) {
      s_l = Hmisc::Lag(sl[, k])
      res = nnet::multinom(s[s_l == 1, i] ~ x[s_l == 1], model = TRUE)


      if(m1 > 2){
        px1 = matrix(stats::fitted(res),
                     ncol = ncol(stats::fitted(res)),
                     nrow = nrow(stats::fitted(res)))
      }else{
        px1 = matrix(cbind(stats::fitted(res), 1-stats::fitted(res)),
                     ncol = ncol(stats::fitted(res))+1,
                     nrow = nrow(stats::fitted(res)))
      }

      for (n in 1:ncol(px1)) {
        if(mov == TRUE){
          ma = pracma::movavg(px1[, n], n = order, type = "s")
          df = data.frame(px = ma, time = seq(1, nrow(px1)))
        }else{
          df = data.frame(px = px1[, n], time = seq(1, nrow(px1)))
        }

        assign(
          paste('p', n, k, paste(i, j, sep = ''), sep = "_"),
          ggplot2::ggplot(df, ggplot2::aes(x = time,
                         y = px)) +
            ggplot2::ylab(label = paste(
              'P(',
              paste('S_', i, 't', sep = ''),
              '=',
              n,
              '|',
              paste('S_', j, 't-1', sep = ''),
              '=',
              k,
              ')'
            )) +
            ggplot2::geom_line(color = 'black') +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.title = ggplot2::element_text(size = 6))
        )
        out <- c(out, paste('p', n, k, paste(i, j, sep = ''), sep = "_"))
      }

    }

    return = plots(out=out)

    ggplot2::ggsave(
      paste(paste('Plot', paste('S', i, j, sep = ''), sep='_'), '.png', sep = ''),
      eval(parse(text= return)),
      width = wdth,
      height = hght,
      units = unts
    )

    out = c()

 }

}


}
