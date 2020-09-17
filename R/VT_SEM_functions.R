perform_SEM_model0=function(dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel,max.cores=1){
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7

  nv=1 #Number of variants
  dat_SEM <- dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel)]
  colnames(dat_SEM)=c("NTm","Tm","NTp","Tp","Yo")

  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=T,values=0,label="F",name="f")
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=T,values=.3,label="EnvPath",name="e")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=T,values=0,label="additive_coefficient",name="delta")
  w <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=T,values=0,label="covAandF",name="w")
  k <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=FALSE,values=0.5,label="varTPs",name="k")
  x1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=.20,label="LatentF1",name="x1")


  # Matrices for variances of phenotypic variances
  sigma <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=T,values=1,label="VarYo",name="sigma")

  # mxAlgebra - nonlinear constraints
  x2 <- mxAlgebra(2*f^2*sigma,name="x2")
  sigma2 <- mxAlgebra(2*delta^2*k + x1 + 2*delta*w +e^2,name="sigma2")
  w2 <- mxAlgebra(2*f/(1-f)*k*delta,name="w2")

  # Equating nonlinear constraints and parameters
  sigmaCon <- mxConstraint(sigma==sigma2,name='sigmaCon')
  xCon <- mxConstraint(x1==x2,name='FCon')
  wCon <-mxConstraint(w==w2,name="wCon")

  # mxAlgebra - relative covariances
  ThetaNTM <- mxAlgebra(.5*f*(w+delta),name="ThetaNTM")
  ThetaTM <- mxAlgebra(.5*f*(w+delta)+k*delta,name="ThetaTM")


  #Put these relative covariances together into MZ relatives and DZ relatives matrices
  CVmat<-    mxAlgebra(rbind(
    cbind(k       ,0           ,0           ,0          ,ThetaNTM   ),
    cbind(0       ,k           ,0           ,0          ,ThetaTM    ),
    cbind(0       ,0           ,k           ,0          ,ThetaNTM   ),
    cbind(0       ,0           ,0           ,k          ,ThetaTM    ),
    cbind(ThetaNTM ,ThetaTM    ,ThetaNTM     ,ThetaTM   ,sigma       )
  ),
  dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")


  # Put in the data
  dat_SEM2 <- mxData( observed=dat_SEM, type="raw" )


  #Objective object
  means <-   mxMatrix(type="Full", nrow=1, ncol=5, free=TRUE,  values= rep(.05,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Yo")), name="expMean")
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Yo"))

  funML <- mxFitFunctionML()
  #params <- list(f,e,delta,w,Vp1,VF1,Vp2,VpCon,CVmat,funML,means,objdat,VFCon,VF2,
  #               w_delta1,w_delta2,w_deltaCon,CvYpYo1,CvYpYo2,CvYpYoCon,CvYmYo1,CvYmYo2,CvYmYoCon,w2,wCon)

  params <- list(f, e , delta, w , k, x1 ,
                 sigma, x2 , sigma2,  w2 ,
                 sigmaCon , xCon, wCon,ThetaNTM , ThetaTM ,
                 CVmat,
                 funML, means, objdat
  )
  init_list=data.frame(f.init=0,e.init=0.3,delta.init=0)

  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  #Run the model
  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  aa=warnings()

  if(sum(str_detect(names(aa),"Mx status RED"))==0&class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k), AFE.Fit$AFEmodel, T)
    VF=mxEval(x1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigmaest=mxEval(sigma,AFE.Fit$AFEmodel,T)
    west=mxEval(w,AFE.Fit$AFEmodel,T)
    h2=VAO/sigmaest
    VFratio=VF/sigmaest
    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=NA
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=NA
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results22=data.frame(method="model0",VAO=VAO,VAL=NA,VF=VF,VE=VE,sigma=sigmaest,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)

  }else{
    results22=data.frame(method="model0",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)


  }
  for(f.init in c(0,0.2,0.4,0.6)){
    for(e.init in c(0.3,0.9)){
      for(delta.init in c(0,0.3,0.7,0.9)){
        init_list=rbind(init_list,data.frame(f.init=f.init,e.init=e.init,delta.init=delta.init))

        f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
        e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
        delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
        params <- list(f, e , delta, w , k, x1 ,
                       sigma, x2 , sigma2,  w2 ,
                       sigmaCon , xCon, wCon,ThetaNTM , ThetaTM ,
                       CVmat,
                       funML, means, objdat
        )

        modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
        AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
        if(class(AFE.Fit)!="try-error"){
          VAO=mxEval(delta^2*2*(k), AFE.Fit$AFEmodel, T)
          VF=mxEval(x1, AFE.Fit$AFEmodel, T)
          VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
          sigmaest=mxEval(sigma,AFE.Fit$AFEmodel,T)
          west=mxEval(w,AFE.Fit$AFEmodel,T)
          h2=VAO/sigmaest
          VFratio=VF/sigmaest
          k_est=mxEval(k,AFE.Fit$AFEmodel,T)
          mu_est=NA
          deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
          aest=NA
          eest=mxEval(e,AFE.Fit$AFEmodel,T)
          f_est=mxEval(f,AFE.Fit$AFEmodel,T)
          ll=summary(AFE.Fit)$Minus2LogLikelihood
          results2=data.frame(method="model0",VAO=VAO,VAL=NA,VF=VF,VE=VE,sigma=sigmaest,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                              deltaest=deltaest,aest=aest,eest=eest,ll=ll)
          results22=rbind(results22,results2)
        }else{
          results2=data.frame(method="model0",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                              deltaest=NA,aest=NA,eest=NA,ll=NA)
          results22=rbind(results22,results2)
        }
      }
    }
  }

  if(sum(is.na(results22$ll))!=nrow(results22)){
    results11=results22[which(results22$ll==min(results22$ll,na.rm=TRUE))[1],]
    f.init=init_list$f.init[ which(results22$ll==min(results22$ll,na.rm=TRUE))[1]]
    e.init=init_list$e.init[ which(results22$ll==min(results22$ll,na.rm=TRUE))[1]]
    delta.init=init_list$delta.init[ which(results22$ll==min(results22$ll,na.rm=TRUE))[1]]
    f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
    e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
    delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
    params <- list(f, e , delta, w , k, x1 ,
                   sigma, x2 , sigma2,  w2 ,
                   sigmaCon , xCon, wCon,ThetaNTM , ThetaTM ,
                   CVmat,
                   funML, means, objdat
    )

    modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
    AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)

  }else{
    results11=data.frame(method="model0",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  results=list(AFE.Fit=AFE.Fit,result_summary=results11)
  return(results)
}

perform_SEM_model1_eq=function(dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel,max.cores=1){
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel)]
  colnames(dat_SEM)=c("NTm","Tm","NTp","Tp","Yo")

  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k")

  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.22,label="F",name="f")
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=1,label="EnvPath",name="e")
  g <-mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.05,label="hap_PRS_cov",name="g")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=sqrt(.5),label="additive_coefficient",name="delta")


  #Constraints - I've checked and all 3 of these are required but no others
  #Matrices for latent variables & nonlinearly constrained estimates
  x1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=.20,label="LatentF1",name="x1")
  w1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.30,label="covAandF",name="w1")
  mu1 <-mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.1,label="copath",name="mu1")
  # mxAlgebra - nonlinear constraints
  x2 <- mxAlgebra(2*f^2*(sigma2+sigma2^2*mu1),name="x2")
  w2 <- mxAlgebra(2*f*Omega2+2*f*sigma2*mu1*Omega2,name="w2")
  mu2 <- mxAlgebra(g/Omega2^2,name="mu2")
  # Equating nonlinear constraints and parameters
  xCon <-mxConstraint(x1==x2,name="xCon")
  wCon <-mxConstraint(w1==w2,name="wCon")
  muCon <- mxConstraint(mu1==mu2,name="muCon")


  # mxAlgebra for implied parameters and relative covariances
  sigma2 <- mxAlgebra(2*delta*Omega2+delta*w1+x1+e^2,name="sigma2")
  Omega2 <- mxAlgebra(2*delta*g+.5*w1+delta*k,name="Omega2")
  ThetaNTM <- mxAlgebra(2*delta*g + f*Omega2*(1+sigma2*mu1),name="ThetaNTM")
  ThetaTM <- mxAlgebra(delta*k + ThetaNTM,name="ThetaTM")


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    cbind(k+g      ,g          ,g        ,g       ,ThetaNTM   ),
    cbind(g        ,k+g        ,g        ,g       ,ThetaTM    ),
    cbind(g        ,g          ,k+g      ,g       ,ThetaNTM   ),
    cbind(g        ,g          ,g        ,k+g     ,ThetaTM    ),
    cbind(ThetaNTM ,ThetaTM    ,ThetaNTM ,ThetaTM ,sigma2     ) ),
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <-   mxMatrix(type="Full", nrow=1, ncol=5, free=TRUE,  values= rep(.05,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Yo")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Yo"))
  params <- list(k,f,e,g,delta,
                 x1,x2,w1,w2,xCon,wCon,mu1,mu2,muCon,
                 sigma2,Omega2,
                 ThetaNTM,ThetaTM,
                 CVmat,MNmat,
                 funML,objdat)


  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  init_list=data.frame(f.init=0.22,e.init=1,g.init=0.05,delta.init=sqrt(.5))

  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  aa=warnings()
  if(sum(str_detect(names(aa),"Mx status RED"))==0&class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
    VF=mxEval(x1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigma=mxEval(sigma2,AFE.Fit$AFEmodel,T)
    west=mxEval(w1,AFE.Fit$AFEmodel,T)
    h2=VAO/sigma
    VFratio=VF/sigma
    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=NA
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results11=data.frame(method="model1_eq",VAO=VAO,VAL=NA,VF=VF,VE=VE,sigma=sigma,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)

  }else{
    results11=data.frame(method="model1_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)


  }

  for(f.init in c(0,0.2,0.4,0.6)){
    for(e.init in c(0.3,0.9)){
      for(g.init in c(0,0.5)){
        for(delta.init in c(0,0.3,0.7,0.9)){
          f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
          e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
          g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
          delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
          init_list=rbind(init_list,data.frame(f.init=f.init,e.init=e.init,g.init=g.init,delta.init=delta.init))
          params <- list(k,f,e,g,delta,
                         x1,x2,w1,w2,xCon,wCon,mu1,mu2,muCon,
                         sigma2,Omega2,
                         ThetaNTM,ThetaTM,
                         CVmat,MNmat,
                         funML,objdat)

          modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
          AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
          if(class(AFE.Fit)!="try-error"){
            VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
            VF=mxEval(x1, AFE.Fit$AFEmodel, T)
            VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
            sigma=mxEval(sigma2,AFE.Fit$AFEmodel,T)
            west=mxEval(w1,AFE.Fit$AFEmodel,T)
            h2=VAO/sigma
            VFratio=VF/sigma
            k_est=mxEval(k,AFE.Fit$AFEmodel,T)
            mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
            deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
            aest=NA
            eest=mxEval(e,AFE.Fit$AFEmodel,T)
            f_est=mxEval(f,AFE.Fit$AFEmodel,T)
            ll=summary(AFE.Fit)$Minus2LogLikelihood

            results1=data.frame(method="model1_eq",VAO=VAO,VAL=NA,VF=VF,VE=VE,sigma=sigma,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                                deltaest=deltaest,aest=aest,eest=eest,ll=ll)
            results11=rbind(results11,results1)
          }else{


            h2true=var(dat$AO)/var(dat$Y)
            results1=data.frame(method="model1_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                                deltaest=NA,aest=NA,eest=NA,ll=NA)
            results11=rbind(results11,results1)
          }
        }
      }
    }
  }


  if(sum(is.na(results11$ll))!=nrow(results11)){
    results11=results11[which(results11$ll==min(results11$ll,na.rm=TRUE))[1],]
    f.init=init_list$f.init[ which(results11$ll==min(results11$ll,na.rm=TRUE))[1]]
    e.init=init_list$e.init[ which(results11$ll==min(results11$ll,na.rm=TRUE))[1]]
    g.init=init_list$g.init[ which(results11$ll==min(results11$ll,na.rm=TRUE))[1]]
    delta.init=init_list$delta.init[ which(results11$ll==min(results11$ll,na.rm=TRUE))[1]]
    f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
    e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
    g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
    delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
    params <- list(k,f,e,g,delta,
                   x1,x2,w1,w2,xCon,wCon,mu1,mu2,muCon,
                   sigma2,Omega2,
                   ThetaNTM,ThetaTM,
                   CVmat,MNmat,
                   funML,objdat)

    modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
    AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)

  }else{
    cat("All iterations failed to be converged\\")
    results11=data.frame(method="model1_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  results=list(AFE.Fit=AFE.Fit,result_summary=results11)

  return(results)
}

perform_SEM_model1_dis=function(dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel,max.cores=1){
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel)]
  colnames(dat_SEM)=c("NTm","Tm","NTp","Tp","Yo")

  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_PRS_t0",name="k") #for var(hap_PRS) if the haplotypic PRS is scaled at equilibrium to have var(hap_PRS) = 1/2
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  #k <- mxAlgebra(.5-g,name="k") #for var(hap_PRS) if the haplotypic PRS is scaled at equilibrium to have var(hap_PRS) = 1/2
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.22,label="F",name="f")
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=1,label="EnvPath",name="e")
  g <-mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.05,label="hap_PRS_cov",name="g")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=sqrt(.5),label="additive_coefficient",name="delta")


  #Constraints - I've checked and all 3 of these are required but no others
  #Matrices for latent variables & nonlinearly constrained estimates
  xo1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=.20,label="LatentFo1",name="xo1")
  xp1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=.20,label="LatentFp1",name="xp1")
  wp1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.30,label="covAandFp",name="wp1")
  wo1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.30,label="covAandFo",name="wo1")
  mu1 <-mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.1,label="copath",name="mu1")
  # mxAlgebra - nonlinear constraints
  xo2 <- mxAlgebra(2*f^2*(sigmap2+sigmap2^2*mu1),name="xo2")
  xp2 <- mxAlgebra(2*f^2*(sigmap2),name="xp2")
  wo2 <- mxAlgebra(2*f*Omega2+2*f*sigmap2*mu1*Omega2,name="wo2")
  wp2 <- mxAlgebra(2*f*Omega2,name="wp2")
  mu2 <- mxAlgebra(g/Omega2^2,name="mu2")
  # Equating nonlinear constraints and parameters
  xoCon <-mxConstraint(xo1==xo2,name="xoCon")
  xpCon <-mxConstraint(xp1==xp2,name="xpCon")
  woCon <-mxConstraint(wo1==wo2,name="woCon")
  wpCon <-mxConstraint(wp1==wp2,name="wCon")
  muCon <- mxConstraint(mu1==mu2,name="muCon")


  # mxAlgebra for implied parameters and relative covariances
  sigmao2 <- mxAlgebra(2*delta^2*k+2*delta^2*g+2*delta*wo1+xo1+e^2,name="sigmao2")
  sigmap2 <- mxAlgebra(2*delta^2*k+2*delta*wp1+xp1+e^2,name="sigmap2")
  Omega2 <- mxAlgebra(.5*wp1+delta*k,name="Omega2")
  ThetaTM05 <- mxAlgebra(delta*k+delta*g+.5*wo1,name="ThetaTM05")
  ThetaNTM05 <- mxAlgebra(delta*g+.5*wo1,name="ThetaNTM05") #this is HALF theta_NTM in the math


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    cbind(k        ,0          ,g        ,g       ,ThetaNTM05   ),
    cbind(0        ,k          ,g        ,g       ,ThetaTM05    ),
    cbind(g        ,g          ,k        ,0       ,ThetaNTM05   ),
    cbind(g        ,g          ,0        ,k       ,ThetaTM05    ),
    cbind(ThetaNTM05 ,ThetaTM05    ,ThetaNTM05 ,ThetaTM05 ,sigmao2     ) ),
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <-   mxMatrix(type="Full", nrow=1, ncol=5, free=TRUE,  values= rep(.05,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Yo")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Yo"))
  params <- list(k ,f ,e , g , delta,
                 xo1, xp1 , wp1 ,wo1 ,mu1 ,
                 xo2 , xp2 ,wo2 ,wp2 ,mu2 ,
                 xoCon, xpCon , woCon , wpCon , muCon ,
                 sigmao2 ,  sigmap2, Omega2,  ThetaTM05, ThetaNTM05,
                 CVmat,MNmat,
                 funML,objdat)


  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  init_list=data.frame(f.init=0.22,e.init=1,g.init=0.05,delta.init=sqrt(.5))

  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  aa=warnings()
  if(sum(str_detect(names(aa),"Mx status RED"))==0&class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+g), AFE.Fit$AFEmodel, T)
    VF=mxEval(xo1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigmao=mxEval(sigmao2,AFE.Fit$AFEmodel,T)
    west=mxEval(wo1,AFE.Fit$AFEmodel,T)
    h2=VAO/sigmao
    VFratio=VF/sigmao
    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=NA
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood

    results66=data.frame(method="model1_dis",VAO=VAO,VAL=NA,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)

  }else{


    h2true=var(dat$AO)/var(dat$Y)
    results66=data.frame(method="model1_dis",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }

  for(f.init in c(0,0.2,0.4,0.6)){
    for(e.init in c(0.3,0.9)){
      for(g.init in c(0,0.5)){
        for(delta.init in c(0,0.3,0.7,0.9)){
          init_list=rbind(init_list,data.frame(f.init=f.init,e.init=e.init,g.init=g.init,delta.init=delta.init))
          f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
          e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
          g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
          delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
          params <- list(k ,f ,e , g , delta,
                         xo1, xp1 , wp1 ,wo1 ,mu1 ,
                         xo2 , xp2 ,wo2 ,wp2 ,mu2 ,
                         xoCon, xpCon , woCon , wpCon , muCon ,
                         sigmao2 ,  sigmap2, Omega2,  ThetaTM05, ThetaNTM05,
                         CVmat,MNmat,
                         funML,objdat)

          modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
          AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
          if(class(AFE.Fit)!="try-error"){
            VAO=mxEval(delta^2*2*(k+g), AFE.Fit$AFEmodel, T)
            VF=mxEval(xo1, AFE.Fit$AFEmodel, T)
            VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
            sigmao=mxEval(sigmao2,AFE.Fit$AFEmodel,T)
            west=mxEval(wo1,AFE.Fit$AFEmodel,T)
            h2=VAO/sigmao
            VFratio=VF/sigmao




            k_est=mxEval(k,AFE.Fit$AFEmodel,T)
            mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
            deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
            aest=NA
            eest=mxEval(e,AFE.Fit$AFEmodel,T)
            f_est=mxEval(f,AFE.Fit$AFEmodel,T)
            ll=summary(AFE.Fit)$Minus2LogLikelihood
            results6=data.frame(method="model1_dis",VAO=VAO,VAL=NA,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                                deltaest=deltaest,aest=aest,eest=eest,ll=ll)
            results66=rbind(results66,results6)


          }else{


            h2true=var(dat$AO)/var(dat$Y)
            results6=data.frame(method="model1_dis",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                                deltaest=NA,aest=NA,eest=NA,ll=NA)
            results66=rbind(results66,results6)
          }
        }
      }
    }
  }

  if(sum(is.na(results66$ll))!=nrow(results66)){
    results11=results66[which(results66$ll==min(results66$ll,na.rm=TRUE))[1],]
    f.init=init_list$f.init[ which(results66$ll==min(results66$ll,na.rm=TRUE))[1]]
    e.init=init_list$e.init[ which(results66$ll==min(results66$ll,na.rm=TRUE))[1]]
    g.init=init_list$g.init[ which(results66$ll==min(results6$ll,na.rm=TRUE))[1]]
    delta.init=init_list$delta.init[ which(results66$ll==min(results66$ll,na.rm=TRUE))[1]]
    f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
    e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
    g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
    delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
    params <- list(k ,f ,e , g , delta,
                   xo1, xp1 , wp1 ,wo1 ,mu1 ,
                   xo2 , xp2 ,wo2 ,wp2 ,mu2 ,
                   xoCon, xpCon , woCon , wpCon , muCon ,
                   sigmao2 ,  sigmap2, Omega2,  ThetaTM05, ThetaNTM05,
                   CVmat,MNmat,
                   funML,objdat)

    modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
    AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)

  }else{
    results11=data.frame(method="model1_dis",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  results=list(AFE.Fit=AFE.Fit,result_summary=results11)
  return(results)
}

perform_SEM_model2_eq=function(dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel,max.cores=1){
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel)]


  #Colnames
  colnames(dat_SEM) <- c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")
  nv=1

  #Scaling factor of the PRS - change depending on how PRS is scaled
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  j <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_latent_hap_PRS_t0",name="j") #for var(hap_PRS_lat)=.5 at t0
  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS_lat)=.5 at t0
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="F",name="f",lbound=-.01)
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.7,label="EnvPath",name="e",lbound=-.01)
  g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="hap_PRS_cov",name="g")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=sqrt(.5),label="obs_coef",name="delta",lbound=-.01)
  a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="latent_coef",name="a",lbound=-.01)


  #Constraints
  #Matrices for latent variables & nonlinearly constrained estimates
  x1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0,label="LatentF1",name="x1")
  w1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covAandF",name="w1")
  mu1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="copath",name="mu1")
  i1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covOL",name="i1")
  #sigma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=var(dat_SEM$Yo),label="sigma",name="sigma1")
  Omega1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Omega",name="Omega1")
  Gamma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Gamma",name="Gamma1")


  # mxAlgebra - nonlinear constraints
  x2 <- mxAlgebra(2*f^2*sigma1*(1 + sigma1*mu1),name="x2")
  w2 <- mxAlgebra(2*f*Omega1*(1 + sigma1*mu1),name="w2")
  v1 <- mxAlgebra((w1*a)/delta,name="v1")
  i2 <- mxAlgebra(Gamma1*mu1*Omega1,name="i2")
  g2 <- mxAlgebra(Omega1^2*mu1,name="g2")
  Omega2 <- mxAlgebra(2*a*i1 + delta*k + 2*delta*g+.5*w1,name="Omega2") #more general: between parent Y and their PRS
  Gamma2 <- mxAlgebra(2*a*h + a*j + 2*delta*i1 + .5*v1,name="Gamma2")
  h <- mxAlgebra((g*a^2)/delta^2,name="h")


  # Equating nonlinear constraints and parameters
  xCon <- mxConstraint(x1==x2,name="xCon")
  wCon <- mxConstraint(w1==w2,name="wCon")
  iCon <- mxConstraint(i1==i2,name="iCon") #not sure if needed - yes, seems to be needed
  gCon <- mxConstraint(g==g2,name="gCon")
  sigmaCon <- mxConstraint(sigma1==sigma2,name="sigmaCon")
  OmegaCon <- mxConstraint(Omega1==Omega2,name="OmegaCon")
  GammaCon <- mxConstraint(Gamma1==Gamma2,name="GammaCon")


  # mxAlgebra for implied parameters and relative covariances
  sigma1 <- mxAlgebra(2*a*Gamma1 + 2*delta*Omega1 + a*v1 + delta*w1 + x1 + e^2,name="sigma1")
  ThetaNTM05 <- mxAlgebra((Omega1-delta*k),name="ThetaNTM05") #this is HALF theta_NTM in the math
  spsPRS <-mxAlgebra(sigma1*mu1*Omega1,name="spsPRS") #THIS HAS BEEN FIXED; more general: between parent Y and spouse's PRS
  PO <- mxAlgebra((a*Gamma1 + delta*Omega1 + f*sigma1)*(1+mu1*sigma1),name="PO")
  spouse <-mxAlgebra(mu1*sigma1^2,name="spouse")


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    #       NTm       Tm          NTp     Tp          Ym         Yp        Yo
    cbind(  k+g        ,g         ,g          ,g        ,Omega1     ,spsPRS    ,ThetaNTM05   ),    #NTm
    cbind(  g          ,k+g       ,g          ,g        ,Omega1     ,spsPRS    ,Omega1        ),    #Tm
    cbind(  g          ,g         ,k+g        ,g        ,spsPRS     ,Omega1    ,ThetaNTM05   ),    #NTp
    cbind(  g          ,g         ,g          ,k+g      ,spsPRS     ,Omega1    ,Omega1     ),    #Tp
    cbind(  Omega1     ,Omega1    ,spsPRS     ,spsPRS   ,sigma1     ,spouse    ,PO         ),    #Ym
    cbind(  spsPRS     ,spsPRS    ,Omega1     ,Omega1   ,spouse     ,sigma1    ,PO         ),    #Yp
    cbind(  ThetaNTM05 ,Omega1    ,ThetaNTM05 ,Omega1   ,PO         ,PO        ,sigma1     ) ),  #Yo
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <- mxMatrix(type="Full", nrow=1, ncol=7, free=TRUE,  values= rep(.01,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYm","meanYp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Ym","Yp","Yo"))

  params <- list(k ,j,  f , e , g ,delta,a,
                 x1 , w1 ,mu1 , i1, sigma1, Omega1 , Gamma1,
                 x2, w2, v1, i2, g2, Omega2 , Gamma2,
                 xCon, wCon, iCon, gCon, OmegaCon, GammaCon,
                 ThetaNTM05, h,
                 spsPRS, PO, spouse,
                 CVmat, MNmat,
                 funML,objdat)


  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  if(class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
    VAL=mxEval(a^2*2*(j+2*h), AFE.Fit$AFEmodel, T)
    VF=mxEval(x1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigma=mxEval(sigma1,AFE.Fit$AFEmodel,T)
    west=mxEval(w1,AFE.Fit$AFEmodel,T)
    h2=(VAO+VAL)/sigma
    VFratio=VF/sigma





    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=mxEval(a,AFE.Fit$AFEmodel,T)
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results44=data.frame(method="model2_eq",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigma,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)
  }else{
    results44=data.frame(method="model2_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }

  init_list=data.frame(f.init=0,e.init=0.7,g.init=0,delta.init=sqrt(.5),a.init=0)
  for(f.init in c(0,0.2,0.4)){
    for(e.init in c(0.3,0.6)){
      for(g.init in c(0,0.5)){
        for(delta.init in c(0,0.3,0.7,0.9)){
          for(a.init in c(0,0.3,0.7,0.9)){
            f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
            e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
            g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
            delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
            a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=a.init,label="latent_additive_coefficient",name="a",lbound=-0.01)
            init_list=rbind(init_list,data.frame(f.init=f.init,e.init=e.init,g.init=g.init,delta.init=delta.init,a.init=a.init))
            params <- list(k ,j,  f , e , g ,delta,a,
                           x1 , w1 ,mu1 , i1, sigma1, Omega1 , Gamma1,
                           x2, w2, v1, i2, g2, Omega2 , Gamma2,
                           xCon, wCon, iCon, gCon, OmegaCon, GammaCon,
                           ThetaNTM05, h,
                           spsPRS, PO, spouse,
                           CVmat, MNmat,
                           funML,objdat)

            modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
            AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
            if(class(AFE.Fit)!="try-error"){
              VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
              VAL=mxEval(a^2*2*(j+2*h), AFE.Fit$AFEmodel, T)
              VF=mxEval(x1, AFE.Fit$AFEmodel, T)
              VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
              sigma=mxEval(sigma1,AFE.Fit$AFEmodel,T)
              west=mxEval(w1,AFE.Fit$AFEmodel,T)
              h2=(VAO+VAL)/sigma
              VFratio=VF/sigma





              k_est=mxEval(k,AFE.Fit$AFEmodel,T)
              mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
              f_est=mxEval(f,AFE.Fit$AFEmodel,T)
              deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
              aest=mxEval(a,AFE.Fit$AFEmodel,T)
              eest=mxEval(e,AFE.Fit$AFEmodel,T)
              ll=summary(AFE.Fit)$Minus2LogLikelihood
              results4=data.frame(method="model2_eq",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigma,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                                  deltaest=deltaest,aest=aest,eest=eest,ll=ll)
              results44=rbind(results44,results4)
            }else{
              results4=data.frame(method="model2_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                                  deltaest=NA,aest=NA,eest=NA,ll=NA)
              results44=rbind(results44,results4)
            }

          }
        }
      }
    }
  }
  if(sum(is.na(results44$ll))!=nrow(results44)){
    results11=results44[which(results44$ll==min(results44$ll,na.rm=TRUE))[1],]
    f.init=init_list$f.init[ which(results44$ll==min(results44$ll,na.rm=TRUE))[1]]
    e.init=init_list$e.init[ which(results44$ll==min(results44$ll,na.rm=TRUE))[1]]
    g.init=init_list$g.init[ which(results44$ll==min(results44$ll,na.rm=TRUE))[1]]
    delta.init=init_list$delta.init[ which(results44$ll==min(results44$ll,na.rm=TRUE))[1]]
    a.init=init_list$a.init[ which(results44$ll==min(results44$ll,na.rm=TRUE))[1]]
    f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
    e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
    g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
    delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
    a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=a.init,label="latent_additive_coefficient",name="a",lbound=-0.01)
                    params <- list(k ,j,  f , e , g ,delta,a,
                                   x1 , w1 ,mu1 , i1, sigma1, Omega1 , Gamma1,
                                   x2, w2, v1, i2, g2, Omega2 , Gamma2,
                                   xCon, wCon, iCon, gCon, OmegaCon, GammaCon,
                                   ThetaNTM05, h,
                                   spsPRS, PO, spouse,
                                   CVmat, MNmat,
                                   funML,objdat)
    AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)

  }else{
    cat("All iterations failed to be converged\\")
    results11=data.frame(method="model2_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                       deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  results=list(AFE.Fit=AFE.Fit,result_summary=results11)

  return(results)
}

perform_SEM_model2_dis=function(dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel,max.cores=1){
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel)]
  colnames(dat_SEM) <- c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")

  #Scaling factor of the PRS - change depending on how PRS is scaled
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  j <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_latent_hap_PRS_t0",name="j") #for var(hap_PRS_lat)=.5 at t0
  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_PRS_t0",name="k") #for var(hap_PRS) if the haplotypic PRS is scaled at equilibrium to have var(hap_PRS) = 1/2
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="F",name="f",lbound=-.01)
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.7,label="EnvPath",name="e",lbound=-.01)
  g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="hap_PRS_cov",name="g")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=sqrt(.5),label="obs_coef",name="delta",lbound=-.01)
  a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="latent_coef",name="a",lbound=-.01)


  #Constraints
  #Matrices for latent variables & nonlinearly constrained estimates
  xo1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0,label="LatentFo1",name="xo1")
  xp1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0,label="LatentFp1",name="xp1")
  wo1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covAandFo",name="wo1")
  wp1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covAandFp",name="wp1")
  mu1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="copath",name="mu1")
  #i1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covOL",name="i1")
  #sigma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=var(dat_SEM$Yo),label="sigma",name="sigma1")
  Omega1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Omega",name="Omega1")
  Gamma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Gamma",name="Gamma1")


  # mxAlgebra - nonlinear constraints
  xo2 <- mxAlgebra(2*f^2*sigmap1*(1 + sigmap1*mu1),name="xo2")
  xp2 <- mxAlgebra(2*f^2*sigmap1,name="xp2")
  wo2 <- mxAlgebra(2*f*Omega1*(1 + sigmap1*mu1),name="wo2")
  wp2 <- mxAlgebra(2*f*Omega1,name="wp2")
  vp <- mxAlgebra(2*f*Gamma1,name="vp")
  vo <- mxAlgebra(2*f*Gamma1+2*f*sigmap1*mu1*Gamma1,name="vo")
  #i2 <- mxAlgebra(Gamma1*mu1*Omega1,name="i2")
  g2 <- mxAlgebra(Omega1^2*mu1,name="g2")
  Omega2 <- mxAlgebra(delta*k +.5*wp1,name="Omega2") #more general: between parent Y and their PRS
  Gamma2 <- mxAlgebra( a*j + .5*vp,name="Gamma2")
  h <- mxAlgebra(Gamma1^2*mu1,name="h")


  # Equating nonlinear constraints and parameters
  xoCon <- mxConstraint(xo1==xo2,name="xoCon")
  xpCon <- mxConstraint(xp1==xp2,name="xpCon")
  woCon <- mxConstraint(wo1==wo2,name="woCon")
  wpCon <- mxConstraint(wp1==wp2,name="wpCon")
  #iCon <- mxConstraint(i1==i2,name="iCon") #not sure if needed - yes, seems to be needed
  gCon <- mxConstraint(g==g2,name="gCon")
  #sigmaCon <- mxConstraint(sigma1==sigma2,name="sigmaCon")
  OmegaCon <- mxConstraint(Omega1==Omega2,name="OmegaCon")
  GammaCon <- mxConstraint(Gamma1==Gamma2,name="GammaCon")


  # mxAlgebra for implied parameters and relative covariances
  sigmao1 <- mxAlgebra(2*a^2*j+2*a^2*h+2*delta^2*k+2*delta^2*g+2*a*vo+2*delta*wo1+xo1+e^2,name="sigmao1")
  sigmap1 <- mxAlgebra(2*a^2*j+2*delta^2*k+2*a*vp+2*delta*wp1+xp1+e^2,name="sigmap1")
  #ThetaTM05 <- mxAlgebra(delta*k+delta*g+f*mu1*sigmap1^2*Omega1+f*Omega1,name="ThetaTM05")
  ThetaTM05 <- mxAlgebra(delta*k+delta*g+a*Gamma1*mu1*Omega1+f*Omega1+f*sigmap1*mu1*Omega1,name="ThetaTM05")
  ThetaNTM05 <- mxAlgebra(delta*g+a*Gamma1*mu1*Omega1+f*Omega1+f*sigmap1*mu1*Omega1,name="ThetaNTM05") #this is HALF theta_NTM in the math
  #ThetaNTM05 <- mxAlgebra(delta*g+f*mu1*sigmap1^2*Omega1+f*Omega1,name="ThetaNTM05") #this is HALF theta_NTM in the math
  spsPRS <-mxAlgebra(sigmap1*mu1*Omega1,name="spsPRS") #THIS HAS BEEN FIXED; more general: between parent Y and spouse's PRS -checked
  PO <- mxAlgebra((a*Gamma1 + delta*Omega1 + f*sigmap1)*(1+mu1*sigmap1),name="PO") #checked
  spouse <-mxAlgebra(mu1*sigmap1^2,name="spouse")# checked


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    #       NTm       Tm          NTp     Tp          Ym         Yp        Yo
    cbind(  k          ,0         ,g          ,g        ,Omega1     ,spsPRS    ,ThetaNTM05   ),    #NTm
    cbind(  0          ,k         ,g          ,g        ,Omega1     ,spsPRS    ,ThetaTM05        ),    #Tm
    cbind(  g          ,g         ,k          ,0        ,spsPRS     ,Omega1    ,ThetaNTM05   ),    #NTp
    cbind(  g          ,g         ,0          ,k        ,spsPRS     ,Omega1    ,ThetaTM05     ),    #Tp
    cbind(  Omega1     ,Omega1    ,spsPRS     ,spsPRS   ,sigmap1    ,spouse    ,PO         ),    #Ym
    cbind(  spsPRS     ,spsPRS    ,Omega1     ,Omega1   ,spouse     ,sigmap1    ,PO         ),    #Yp
    cbind(  ThetaNTM05 ,ThetaTM05    ,ThetaNTM05 ,ThetaTM05  ,PO         ,PO        ,sigmao1     ) ),  #Yo
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <- mxMatrix(type="Full", nrow=1, ncol=7, free=TRUE,  values= rep(.01,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYm","meanYp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Ym","Yp","Yo"))

  params <- list( j , k , f , e ,g ,delta , a,
                  xo1 ,xp1 ,wo1 ,wp1,  mu1,
                  Omega1 , Gamma1 ,
                  xo2, xp2 , wo2, wp2, vp , vo ,
                  g2 ,Omega2,Gamma2 ,h,
                  xoCon ,xpCon, woCon, wpCon ,gCon , OmegaCon ,GammaCon,
                  sigmao1 , sigmap1,ThetaNTM05 , ThetaTM05,spsPRS , PO ,  spouse,
                  CVmat, MNmat,
                  funML,objdat)

  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  if(class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+g), AFE.Fit$AFEmodel, T)
    VAL=mxEval(a^2*2*(j+h), AFE.Fit$AFEmodel, T)
    VF=mxEval(xo1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigmao=mxEval(sigmao1,AFE.Fit$AFEmodel,T)
    sigmpp=mxEval(sigmap1,AFE.Fit$AFEmodel,T)
    west=mxEval(wo1,AFE.Fit$AFEmodel,T)
    h2=(VAO+VAL)/sigmao
    VFratio=VF/sigmao
    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=mxEval(a,AFE.Fit$AFEmodel,T)
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results55=data.frame(method="model2_dis",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)
  }else{
    results55=data.frame(method="model2_dis",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }

  init_list=data.frame(f.init=0,e.init=0.7,g.init=0,delta.init=sqrt(.5),a.init=0)
  for(f.init in c(0,0.2,0.4,0.6)){
    for(e.init in c(0.3,0.9)){
      for(g.init in c(0,0.5)){
        for(delta.init in c(0,0.3,0.7,0.9)){
          for(a.init in c(0,0.3,0.7,0.9)){
            f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
            e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
            g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
            delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
            a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=a.init,label="latent_additive_coefficient",name="a",lbound=-0.01)
            init_list2=data.frame(f.init=f.init,e.init=e.init,g.init=g.init,delta.init=delta.init,a.init=a.init)
            init_list=rbind(init_list,init_list2)
            params <- list( j , k , f , e ,g ,delta , a,
                            xo1 ,xp1 ,wo1 ,wp1,  mu1,
                            Omega1 , Gamma1 ,
                            xo2, xp2 , wo2, wp2, vp , vo ,
                            g2 ,Omega2,Gamma2 ,h,
                            xoCon ,xpCon, woCon, wpCon ,gCon , OmegaCon ,GammaCon,
                            sigmao1 , sigmap1,ThetaNTM05 , ThetaTM05,spsPRS , PO ,  spouse,
                            CVmat, MNmat,
                            funML,objdat)

            modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
            AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
            if(class(AFE.Fit)!="try-error"){
              VAO=mxEval(delta^2*2*(k+g), AFE.Fit$AFEmodel, T)
              VAL=mxEval(a^2*2*(j+h), AFE.Fit$AFEmodel, T)
              VF=mxEval(xo1, AFE.Fit$AFEmodel, T)
              VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
              sigmao=mxEval(sigmao1,AFE.Fit$AFEmodel,T)
              sigmpp=mxEval(sigmap1,AFE.Fit$AFEmodel,T)
              west=mxEval(wo1,AFE.Fit$AFEmodel,T)
              h2=(VAO+VAL)/sigmao
              VFratio=VF/sigmao
              k_est=mxEval(k,AFE.Fit$AFEmodel,T)
              mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
              f_est=mxEval(f,AFE.Fit$AFEmodel,T)
              deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
              aest=mxEval(a,AFE.Fit$AFEmodel,T)
              eest=mxEval(e,AFE.Fit$AFEmodel,T)
              ll=summary(AFE.Fit)$Minus2LogLikelihood
              results5=data.frame(method="model2_dis",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                                  deltaest=deltaest,aest=aest,eest=eest,ll=ll)
              results55=rbind(results55,results5)
            }else{
              results5=data.frame(method="model2_dis",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                                  deltaest=NA,aest=NA,eest=NA,ll=NA)
              results55=rbind(results55,results5)
            }
          }
        }
      }
    }
  }
  if(sum(is.na(results55$ll))!=nrow(results55)){
    results11=results55[which(results55$ll==min(results55$ll,na.rm=TRUE))[1],]
    f.init=init_list$f.init[ which(results55$ll==min(results55$ll,na.rm=TRUE))[1]]
    e.init=init_list$e.init[ which(results55$ll==min(results55$ll,na.rm=TRUE))[1]]
    g.init=init_list$g.init[ which(results55$ll==min(results55$ll,na.rm=TRUE))[1]]
    delta.init=init_list$delta.init[ which(results55$ll==min(results55$ll,na.rm=TRUE))[1]]
    a.init=init_list$a.init[ which(results55$ll==min(results55$ll,na.rm=TRUE))[1]]
    f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
    e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
    g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
    delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
    a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=a.init,label="latent_additive_coefficient",name="a",lbound=-0.01)

    params <- list( j , k , f , e ,g ,delta , a,
                    xo1 ,xp1 ,wo1 ,wp1,  mu1,
                    Omega1 , Gamma1 ,
                    xo2, xp2 , wo2, wp2, vp , vo ,
                    g2 ,Omega2,Gamma2 ,h,
                    xoCon ,xpCon, woCon, wpCon ,gCon , OmegaCon ,GammaCon,
                    sigmao1 , sigmap1,ThetaNTM05 , ThetaTM05,spsPRS , PO ,  spouse,
                    CVmat, MNmat,
                    funML,objdat)
    dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
    modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
    AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)

  }else{
    cat("All iterations failed to be converged\\")
    results11=data.frame(method="model2_dis",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                       deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  results=list(AFE.Fit=AFE.Fit,result_summary=results11)

  return(results)
}

perform_SEM_model2_eq_NP=function(dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel,herit0,prop.latent,max.cores=1){
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel)]
  colnames(dat_SEM) <- c("NTm","Tm","NTp","Tp","Yo")

  #Scaling factor of the PRS - change depending on how PRS is scaled
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  j <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_latent_hap_PRS_t0",name="j") #for var(hap_PRS_lat)=.5 at t0
  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_PRS_t0",name="k") #for var(hap_PRS) if the haplotypic PRS is scaled at equilibrium to have var(hap_PRS) = 1/2
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="F",name="f",lbound=-.01)
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.7,label="EnvPath",name="e",lbound=-.01)
  g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="hap_PRS_cov",name="g")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=sqrt(.5),label="obs_coef",name="delta",lbound=-.01)
  a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=FALSE,values=sqrt(herit0*prop.latent),label="latent_coef",name="a")


  #Constraints
  #Matrices for latent variables & nonlinearly constrained estimates
  xo1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0,label="LatentFo1",name="xo1")
  xp1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0,label="LatentFp1",name="xp1")
  wo1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covAandFo",name="wo1")
  wp1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covAandFp",name="wp1")
  mu1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="copath",name="mu1")
  #i1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covOL",name="i1")
  #sigma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=var(dat_SEM$Yo),label="sigma",name="sigma1")
  Omega1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Omega",name="Omega1")
  Gamma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Gamma",name="Gamma1")


  # mxAlgebra - nonlinear constraints
  xo2 <- mxAlgebra(2*f^2*sigmap1*(1 + sigmap1*mu1),name="xo2")
  xp2 <- mxAlgebra(2*f^2*sigmap1,name="xp2")
  wo2 <- mxAlgebra(2*f*Omega1*(1 + sigmap1*mu1),name="wo2")
  wp2 <- mxAlgebra(2*f*Omega1,name="wp2")
  vp <- mxAlgebra(2*f*Gamma1,name="vp")
  vo <- mxAlgebra(2*f*Gamma1+2*f*sigmap1*mu1*Gamma1,name="vo")
  #i2 <- mxAlgebra(Gamma1*mu1*Omega1,name="i2")
  g2 <- mxAlgebra(Omega1^2*mu1,name="g2")
  Omega2 <- mxAlgebra(delta*k +.5*wp1,name="Omega2") #more general: between parent Y and their PRS
  Gamma2 <- mxAlgebra( a*j + .5*vp,name="Gamma2")
  h <- mxAlgebra(Gamma1^2*mu1,name="h")


  # Equating nonlinear constraints and parameters
  xoCon <- mxConstraint(xo1==xo2,name="xoCon")
  xpCon <- mxConstraint(xp1==xp2,name="xpCon")
  woCon <- mxConstraint(wo1==wo2,name="woCon")
  wpCon <- mxConstraint(wp1==wp2,name="wpCon")
  #iCon <- mxConstraint(i1==i2,name="iCon") #not sure if needed - yes, seems to be needed
  gCon <- mxConstraint(g==g2,name="gCon")
  #sigmaCon <- mxConstraint(sigma1==sigma2,name="sigmaCon")
  OmegaCon <- mxConstraint(Omega1==Omega2,name="OmegaCon")
  GammaCon <- mxConstraint(Gamma1==Gamma2,name="GammaCon")


  # mxAlgebra for implied parameters and relative covariances
  sigmao1 <- mxAlgebra(2*a^2*j+2*a^2*h+2*delta^2*k+2*delta^2*g+2*a*vo+2*delta*wo1+xo1+e^2,name="sigmao1")
  sigmap1 <- mxAlgebra(2*a^2*j+2*delta^2*k+2*a*vp+2*delta*wp1+xp1+e^2,name="sigmap1")
  #ThetaTM05 <- mxAlgebra(delta*k+delta*g+f*mu1*sigmap1^2*Omega1+f*Omega1,name="ThetaTM05")
  ThetaTM05 <- mxAlgebra(delta*k+delta*g+a*Gamma1*mu1*Omega1+f*Omega1+f*sigmap1*mu1*Omega1,name="ThetaTM05")
  ThetaNTM05 <- mxAlgebra(delta*g+a*Gamma1*mu1*Omega1+f*Omega1+f*sigmap1*mu1*Omega1,name="ThetaNTM05") #this is HALF theta_NTM in the math
  #ThetaNTM05 <- mxAlgebra(delta*g+f*mu1*sigmap1^2*Omega1+f*Omega1,name="ThetaNTM05") #this is HALF theta_NTM in the math
  #spsPRS <-mxAlgebra(sigmap1*mu1*Omega1,name="spsPRS") #THIS HAS BEEN FIXED; more general: between parent Y and spouse's PRS -checked
  #PO <- mxAlgebra((a*Gamma1 + delta*Omega1 + f*sigmap1)*(1+mu1*sigmap1),name="PO") #checked
  #spouse <-mxAlgebra(mu1*sigmap1^2,name="spouse")# checked


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    #       NTm       Tm          NTp     Tp          Yo
    cbind(  k          ,0         ,g          ,g       ,ThetaNTM05   ),    #NTm
    cbind(  0          ,k         ,g          ,g       ,ThetaTM05        ),    #Tm
    cbind(  g          ,g         ,k          ,0        ,ThetaNTM05   ),    #NTp
    cbind(  g          ,g         ,0          ,k        ,ThetaTM05     ),    #Tp
    cbind(  ThetaNTM05 ,ThetaTM05    ,ThetaNTM05 ,ThetaTM05 ,sigmao1     ) ),  #Yo
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <- mxMatrix(type="Full", nrow=1, ncol=5, free=TRUE,  values= rep(.01,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Yo")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Yo"))

  params <- list( j , k , f , e ,g ,delta , a,
                  xo1 ,xp1 ,wo1 ,wp1,  mu1,
                  Omega1 , Gamma1 ,
                  xo2, xp2 , wo2, wp2, vp , vo ,
                  g2 ,Omega2,Gamma2 ,h,
                  xoCon ,xpCon, woCon, wpCon ,gCon , OmegaCon ,GammaCon,
                  sigmao1 , sigmap1,ThetaNTM05 , ThetaTM05 ,
                  CVmat, MNmat,
                  funML,objdat)

  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  if(class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+g), AFE.Fit$AFEmodel, T)
    VAL=mxEval(a^2*2*(j+h), AFE.Fit$AFEmodel, T)
    VF=mxEval(xo1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigmao=mxEval(sigmao1,AFE.Fit$AFEmodel,T)
    sigmpp=mxEval(sigmap1,AFE.Fit$AFEmodel,T)
    west=mxEval(wo1,AFE.Fit$AFEmodel,T)
    h2=(VAO+VAL)/sigmao
    VFratio=VF/sigmao




    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)

    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=mxEval(a,AFE.Fit$AFEmodel,T)
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results88=data.frame(method="model2_eq_NP",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)
  }else{
    results88=data.frame(method="model2_eq_NP",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }

  init_list=data.frame(f.init=0.22,e.init=2,g.init=0.005,delta.init=.22)

  for(f.init in c(0,0.2,0.4,0.6)){
    for(e.init in c(0.3,0.9)){
      for(g.init in c(0,0.5)){
        for(delta.init in c(0,0.3,0.7,0.9)){
          init_list=rbind(init_list,data.frame(f.init=f.init,e.init=e.init,g.init=g.init,delta.init=delta.init))
          f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
          e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
          g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
          delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)

          params <- list( j , k , f , e ,g ,delta , a,
                          xo1 ,xp1 ,wo1 ,wp1,  mu1,
                          Omega1 , Gamma1 ,
                          xo2, xp2 , wo2, wp2, vp , vo ,
                          g2 ,Omega2,Gamma2 ,h,
                          xoCon ,xpCon, woCon, wpCon ,gCon , OmegaCon ,GammaCon,
                          sigmao1 , sigmap1,ThetaNTM05 , ThetaTM05 ,
                          CVmat, MNmat,
                          funML,objdat)

          modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
          AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
          if(class(AFE.Fit)!="try-error"){
            VAO=mxEval(delta^2*2*(k+g), AFE.Fit$AFEmodel, T)
            VAL=mxEval(a^2*2*(j+h), AFE.Fit$AFEmodel, T)
            VF=mxEval(xo1, AFE.Fit$AFEmodel, T)
            VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
            sigmao=mxEval(sigmao1,AFE.Fit$AFEmodel,T)
            sigmpp=mxEval(sigmap1,AFE.Fit$AFEmodel,T)
            west=mxEval(wo1,AFE.Fit$AFEmodel,T)
            h2=(VAO+VAL)/sigmao
            VFratio=VF/sigmao





            k_est=mxEval(k,AFE.Fit$AFEmodel,T)
            mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
            f_est=mxEval(f,AFE.Fit$AFEmodel,T)
            deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
            aest=mxEval(a,AFE.Fit$AFEmodel,T)
            eest=mxEval(e,AFE.Fit$AFEmodel,T)
            ll=summary(AFE.Fit)$Minus2LogLikelihood
            results8=data.frame(method="model2_eq_NP",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                                deltaest=deltaest,aest=aest,eest=eest,ll=ll)
            results88=rbind(results88,results8)
          }else{
            results8=data.frame(method="model2_eq_NP",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                                deltaest=NA,aest=NA,eest=NA,ll=NA)
            results88=rbind(results88,results8)
          }
        }
      }
    }
  }

  if(sum(is.na(results88$ll))!=nrow(results88)){
    results11=results88[which(results88$ll==min(results88$ll,na.rm=TRUE))[1],]
    f.init=init_list$f.init[ which(results88$ll==min(results88$ll,na.rm=TRUE))[1]]
    e.init=init_list$e.init[ which(results88$ll==min(results88$ll,na.rm=TRUE))[1]]
    g.init=init_list$g.init[ which(results88$ll==min(results88$ll,na.rm=TRUE))[1]]
    delta.init=init_list$delta.init[ which(results88$ll==min(results88$ll,na.rm=TRUE))[1]]
    f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
    e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
    g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
    delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)

    params <- list( j , k , f , e ,g ,delta , a,
                    xo1 ,xp1 ,wo1 ,wp1,  mu1,
                    Omega1 , Gamma1 ,
                    xo2, xp2 , wo2, wp2, vp , vo ,
                    g2 ,Omega2,Gamma2 ,h,
                    xoCon ,xpCon, woCon, wpCon ,gCon , OmegaCon ,GammaCon,
                    sigmao1 , sigmap1,ThetaNTM05 , ThetaTM05 ,
                    CVmat, MNmat,
                    funML,objdat)
    modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
    AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)

  }else{
    results11=data.frame(method="model2_eq_NP",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  results=list(AFE.Fit=AFE.Fit,result_summary=results11)
  return(results)
}

perform_SEM_model2_dis_NP=function(dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel,herit0,prop.latent,max.cores=1){
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Yolabel)]
  colnames(dat_SEM) <- c("NTm","Tm","NTp","Tp","Yo")

  #Scaling factor of the PRS - change depending on how PRS is scaled
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  j <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_latent_hap_PRS_t0",name="j") #for var(hap_PRS_lat)=.5 at t0
  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS_lat)=.5 at t0
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="F",name="f",lbound=-.01,ubound=1)
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=2,label="EnvPath",name="e",lbound=-.01,ubound=3)
  g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.005,label="hap_PRS_cov",name="g")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.22,label="obs_coef",name="delta",lbound=-.01,ubound=.8)
  a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=FALSE,values=sqrt(herit0*prop.latent),label="latent_coef",name="a")



  #Constraints
  #Matrices for latent variables & nonlinearly constrained estimates
  x1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0,label="LatentF1",name="x1")
  w1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covAandF",name="w1")
  v1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covLatandF",name="v1")
  mu1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.092,label="copath",name="mu1")
  i1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="covOL",name="i1")
  Omega1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Omega",name="Omega1")
  Gamma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0,label="Gamma",name="Gamma1")


  # mxAlgebra - nonlinear constraints
  x2 <- mxAlgebra(2*f^2*sigma1*(1 + sigma1*mu1),name="x2")
  w2 <- mxAlgebra(2*f*Omega1*(1 + sigma1*mu1),name="w2")
  v2 <- mxAlgebra((w1*a)/delta,name="v2")
  i2 <- mxAlgebra(Gamma1*mu1*Omega1,name="i2")
  g2 <- mxAlgebra(Omega1^2*mu1,name="g2")
  Omega2 <- mxAlgebra(2*a*i1 + delta*k + 2*delta*g+.5*w1,name="Omega2") #more general: between parent Y and their PRS
  Gamma2 <- mxAlgebra(2*a*h + a*j + 2*delta*i1 + .5*v1,name="Gamma2")
  h <- mxAlgebra((g*a^2)/delta^2,name="h")
  #mu2 <- mxAlgebra(g/Omega2^2,name="mu2")


  # Equating nonlinear constraints and parameters
  xCon <- mxConstraint(x1==x2,name="xCon")
  wCon <- mxConstraint(w1==w2,name="wCon")
  vCon <- mxConstraint(v1==v2,name="vCon")
  iCon <- mxConstraint(i1==i2,name="iCon") #not sure if needed - yes, seems to be needed
  gCon <- mxConstraint(g==g2,name="gCon")
  OmegaCon <- mxConstraint(Omega1==Omega2,name="OmegaCon")
  GammaCon <- mxConstraint(Gamma1==Gamma2,name="GammaCon")
  #muCon <- mxConstraint(mu1==mu2,name="muCon")


  # mxAlgebra for implied parameters and relative covariances
  sigma1 <- mxAlgebra(2*a*Gamma1 + 2*delta*Omega1 + a*v1 + delta*w1 + x1 + e^2,name="sigma1")
  ThetaNTM05 <- mxAlgebra((Omega1-delta*k),name="ThetaNTM05") #this is HALF theta_NTM in the math
  #spsPRS <-mxAlgebra(sigma1*mu1*Omega1,name="spsPRS") #THIS HAS BEEN FIXED; more general: between parent Y and spouse's PRS
  #PO <- mxAlgebra((a*Gamma1 + delta*Omega1 + f*sigma1)*(1+mu1*sigma1),name="PO")
  #spouse <-mxAlgebra(mu1*sigma1^2,name="spouse")


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    #       NTm       Tm          NTp     Tp              Yo
    cbind(  k+g        ,g         ,g          ,g        ,ThetaNTM05   ),    #NTm
    cbind(  g          ,k+g       ,g          ,g        ,Omega1        ),    #Tm
    cbind(  g          ,g         ,k+g        ,g        ,ThetaNTM05   ),    #NTp
    cbind(  g          ,g         ,g          ,k+g      ,Omega1     ),    #Tp
    cbind(  ThetaNTM05 ,Omega1    ,ThetaNTM05 ,Omega1   ,sigma1     ) ),  #Yo
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <- mxMatrix(type="Full", nrow=1, ncol=5, free=TRUE,  values= rep(.01,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Yo")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Yo"))

  params <- list(k ,j,  f , e , g ,delta,a,
                 x1 , w1 ,v1, mu1, i1, sigma1, Omega1 , Gamma1,
                 x2, w2, v2, i2, g2, Omega2 , Gamma2,
                 xCon, wCon, vCon, iCon, gCon, OmegaCon, GammaCon,
                 ThetaNTM05, h,
                 CVmat, MNmat,
                 funML,objdat)


  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)

  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  aa=warnings()
  if(sum(str_detect(names(aa),"Mx status RED"))==0&class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
    VAL=mxEval(a^2*2*(j+2*h), AFE.Fit$AFEmodel, T)
    VF=mxEval(x1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigmao=mxEval(sigma1,AFE.Fit$AFEmodel,T)
    west=mxEval(w1,AFE.Fit$AFEmodel,T)
    h2=VAO/sigmao
    VFratio=VF/sigmao



    h2true=var(dat$AOy)/var(dat$Y)

    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=1
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results77=data.frame(method="model2_dis_NP",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)

  }else{


    h2true=var(dat$AO)/var(dat$Y)
    results77=data.frame(method="model2_dis_NP",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }

  init_list=data.frame(f.init=0.22,e.init=2,g.init=0.005,delta.init=.22)
  for(f.init in c(0,0.2,0.4,0.6)){
    for(e.init in c(0.3,0.9)){
      for(g.init in c(0,0.5)){
        for(delta.init in c(0,0.3,0.7,0.9)){
          init_list=rbind(init_list,data.frame(f.init=f.init,e.init=e.init,g.init=g.init,delta.init=delta.init))
          f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
          e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
          g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
          delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)
          params <- list(k ,j,  f , e , g ,delta,a,
                         x1 , w1 ,v1, mu1, i1, sigma1, Omega1 , Gamma1,
                         x2, w2, v2, i2, g2, Omega2 , Gamma2,
                         xCon, wCon, vCon, iCon, gCon, OmegaCon, GammaCon,
                         ThetaNTM05, h,
                         CVmat, MNmat,
                         funML,objdat)


          modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
          AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
          if(class(AFE.Fit)!="try-error"){
            VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
            VAL=mxEval(a^2*2*(j+2*h), AFE.Fit$AFEmodel, T)
            VF=mxEval(x1, AFE.Fit$AFEmodel, T)
            VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
            sigmao=mxEval(sigma1,AFE.Fit$AFEmodel,T)
            west=mxEval(w1,AFE.Fit$AFEmodel,T)
            h2=VAO/sigmao
            VFratio=VF/sigmao



            h2true=var(dat$AOy)/var(dat$Y)

            k_est=mxEval(k,AFE.Fit$AFEmodel,T)
            mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
            deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
            aest=mxEval(a, AFE.Fit$AFEmodel, T)
            eest=mxEval(e,AFE.Fit$AFEmodel,T)
            f_est=mxEval(f,AFE.Fit$AFEmodel,T)
            ll=summary(AFE.Fit)$Minus2LogLikelihood
            results7=data.frame(method="model2_dis_NP",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigmao,west=west,h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                                deltaest=deltaest,aest=aest,eest=eest,ll=ll)
            results77=rbind(results77,results7)
          }else{
            results7=data.frame(method="model2_dis_NP",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                                deltaest=NA,aest=NA,eest=NA,ll=NA)
            results77=rbind(results77,results7)
          }
        }
      }
    }
  }

  if(sum(is.na(results77$ll))!=nrow(results77)){
    results11=results77[which(results77$ll==min(results77$ll,na.rm=TRUE))[1],]
    f.init=init_list$f.init[ which(results77$ll==min(results77$ll,na.rm=TRUE))[1]]
    e.init=init_list$e.init[ which(results77$ll==min(results77$ll,na.rm=TRUE))[1]]
    g.init=init_list$g.init[ which(results77$ll==min(results77$ll,na.rm=TRUE))[1]]
    delta.init=init_list$delta.init[ which(results77$ll==min(results77$ll,na.rm=TRUE))[1]]
    f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=f.init,label="F",name="f",lbound=-0.01)
    e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=e.init,label="EnvPath",name="e",lbound=-0.01)
    g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=g.init,label="hap_PRS_cov",name="g")
    delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=delta.init,label="additive_coefficient",name="delta",lbound=-0.01)

    modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
    AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)

  }else{
    cat("All iterations failed to be converged\\")
    results11=data.frame(method="model2_dis_NP",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                       deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  results=list(AFE.Fit=AFE.Fit,result_summary=results11)
  return(results)
}

