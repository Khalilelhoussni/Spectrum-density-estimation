wr_vec=function(alpha.vec, p0,yt, filter.number, family="DaubExPhase", bc="periodic"){
  J=log2(p0)
  n=length(alpha.vec)
  #yta=wd(rep(0,n),filter.number=filter.number, family=family, bc=bc)
  yt=putC(yt,J,alpha.vec[1:p0])
  yt$D[1:(n-p0)]=alpha.vec[-c(1:p0)]
  out=wr(yt, start.level=J)
  return(out)
}
