bootstrap=function(model1,model2){
  B=200
  n=nrow(model2$mu.x)
  #mu.link=model1$mu.link
  #phi.link=model1$sigma.link
  mu.formula.1=as.formula(paste(
    "yb~",model1$mu.formula[3]))
  mu.formula.2=as.formula(paste(
    "yb~",model2$mu.formula[3]))
  phi.formula.1=model1$sigma.formula
  phi.formula.2=model2$sigma.formula
  LRb=NULL
  for(i in 1:B){
    yb=rBP(n,model2$mu.fv,model2$sigma.fv)
    ajuste1=tryCatch(gamlss(mu.formula.1,
                            sigma.formula = phi.formula.1,
                            family=BP(mu.link="log",
                                      sigma.link="log"),
                            control=gamlss.control(trace=F)),
                     error=function(e){NA})
    ajuste2=tryCatch(gamlss(mu.formula.2,
                            sigma.formula = phi.formula.2,
                            family=BP(mu.link="log",
                                      sigma.link="log"),
                            control=gamlss.control(trace=F)),
                     error=function(e){NA})
    LRb=c(LRb,tryCatch(abs(
      ajuste1$G.deviance-ajuste2$G.deviance),
                       error=function(e){NA}))
    #print(i/B)
  }
  LRbbar=mean(LRb,na.rm=T)
  return(LRbbar)
}
