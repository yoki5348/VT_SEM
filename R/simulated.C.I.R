c.i.compute_VF=function(n,r2PGS,VF,percentile){
  logn=log10(n)
  sqrtpgs=sqrt(r2PGS)
  se.comp=t(c(0.4661169,7.864836E-7,-7.549848E-2,2.792433E-1,-4.352591E-1))%*%c(1,n,logn,r2PGS,sqrtpgs)
  result=data.frame(Lower.C.I=VF-se.comp*qnorm(1-(1-percentile)/2),Upper.C.I=VF+se.comp*qnorm(1-(1-percentile)/2))
  print(result)
}

c.i.compute_r2PGS=function(n,r2PGS,percentile){
  logn=log10(n)
  sqrtpgs=sqrt(r2PGS)
  se.comp=t(c(5.987E-2,1.801E-7,-1.735E-2,-3.917E-2,5.403E-2))%*%c(1,n,logn,r2PGS,sqrtpgs)
  result=data.frame(Lower.C.I=r2PGS-se.comp*qnorm(1-(1-percentile)/2),Upper.C.I=r2PGS+se.comp*qnorm(1-(1-percentile)/2))
  print(result)
}

c.i.compute_r2LGS=function(n,r2PGS,r2LGS,percentile){
  logn=log10(n)
  sqrtpgs=sqrt(r2PGS)
  se.comp=t(c(5.546E-1,8.066E-7,-7.796E-2,4.865E-1,-7.081E-1))%*%c(1,n,logn,r2PGS,sqrtpgs)
  result=data.frame(Lower.C.I=r2LGS-se.comp*qnorm(1-(1-percentile)/2),Upper.C.I=r2LGS+se.comp*qnorm(1-(1-percentile)/2))
  print(result)
}

c.i.compute_w=function(n,r2PGS,w,percentile){
  logn=log10(n)
  sqrtpgs=sqrt(r2PGS)
  se.comp=t(c(1.171E-1,2.661E-7,-2.616E-2,-3.205E-2,1.197E-2))%*%c(1,n,logn,r2PGS,sqrtpgs)
  result=data.frame(Lower.C.I=w-se.comp*qnorm(1-(1-percentile)/2),Upper.C.I=w+se.comp*qnorm(1-(1-percentile)/2))
  print(result)
}
