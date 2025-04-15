scp <- function(cal, y, test, sig = 0.05) {

  f.score <- function(y, x) abs(y - x)/x
  
  scores <- f.score(y, cal)
  
  n_cal <- length(scores)
  
  q_cal <- quantile(scores, probs = ceiling((1-sig) * (n_cal+1))/n_cal)
  
  ic <- data.frame(test, down = test - q_cal * test, up = test + q_cal * test)
    
  ic$amp <- ic$up - ic$down
  
  return(list(quantile = q_cal, ic_test = ic))
  
}