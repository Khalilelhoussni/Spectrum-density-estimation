wd_vec=function(y, p0, filter.number=1, family="DaubExPhase", bc="periodic"){
  n=length(y)
  J=log2(p0)
  yt=wd(y,filter.number=filter.number, family=family, bc=bc)
  yt.vec=rep(NA,n)
  yt.vec[-c(1:p0)]=yt$D[1:(n-p0)]
  yt.vec[1:p0]=accessC(yt,J)
  return(yt.vec)
}
