ek=function(model){
  mu=model$mu.fv
  phi=model$sigma.fv
  K_1=vcov(model)
  X=model$mu.x
  mu.link=model$mu.link
  phi.link=model$sigma.link
  if(mu.link=="identity"){
    a1=1; a2<-a3<-0
  }
  if(phi.link=="identity"){
    b1=1; b2<-b3<-0
  }
  if(mu.link=="log"){
    a1<-a2<-a3<-mu
  }
  if(phi.link=="log"){
    b1<-b2<-b3<-phi
  }
  if(mu.link=="inverse"){
    a1=-mu^2; a2=2*mu^3; a3=-6*mu^4
  }
  if(phi.link=="inverse"){
    b1=-phi^2; b2=2*phi^3; b3=-6*phi^4
  }
  if(mu.link=="sqrt"){
    a1=2*X%*%model$mu.coef; a2=2; a3 = 0
  }
  if(phi.link=="sqrt"){
    b1=2*model$sigma.coef; b2=2; b3 = 0
  }
  p=ncol(X)
  k=p+1
  ####  Definindo algumas funções ####
  manyzeros=function(r){
    sum(as.numeric(r==0))
  }
  alpha=mu*(1+phi)
  beta=phi+2
  d=function(k,m){
    (1+phi)^k*(psigamma(alpha+beta,m-1)-
                 psigamma(alpha,m-1))
  }
  f=function(k,m){
    (1+phi)^k*psigamma(alpha+beta,m-1)
  }
  j=function(m){
    (1+mu)^m*psigamma(alpha+beta,m-1)-
      mu^m*psigamma(alpha,m-1)-
      psigamma(beta,m-1)
  }
  
 
  
  #### kappas ####
  k4 = function(r,s,t,u){# k_rstu e afins
    vetor_theta = c(r,s,t,u) # vetor de índices do theta 0=phi
    ordem_indices = order(vetor_theta,decreasing=T)
    nphi=manyzeros(vetor_theta) # contando os phis
    
    if(nphi==0){ #k_rstu
      kappa=sum((
        d(4,4)*a1^4+6*d(3,3)*a1^2*a2+d(2,2)*(3*a2^2+4*a1*a3)
      )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]*X[,vetor_theta[ordem_indices[4]]]
      )
    }
    if(nphi==1){ #k_rstphi
      kappa=sum((
        d(3,4)*mu*a1^3+
          3*d(2,3)*(mu*a1*a2+a1^3)+
          d(1,2)*(mu*a3+6*a1*a2)+
          f(3,4)*a1^3+
          3*f(2,3)*a1*a2+
          f(1,2)*a3
      )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
      )
    }
    if(nphi==2){ #k_rsphiphi
      kappa=sum((
        d(2,4)*mu^2*a1^2*b1^2+
          d(2,3)*mu*a1^2*b2+
          d(1,3)*mu*b1^2*(4*a1^2+mu*a2)+
          d(1,2)*b2*(mu*a2+2*a1^2)+
          2*d(0,2)*b1^2*(mu*a2+a1^2)+
          f(2,4)*(1+2*mu)*a1^2*b1^2+
          f(2,3)*a1^2*b2+
          f(1,3)*b1^2*(4*a1^2+(1+2*mu)*a2)+
          f(1,2)*a2*b2+
          2*f(0,2)*a2*b1^2
      )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
      )
    }
    if(nphi==3){ #k_rphiphiphi
      kappa=sum((
        d(1,4)*b1^3*mu^3+
          3*d(1,3)*b1*b2*mu^2+
          d(1,2)*b3*mu+
          3*d(0,3)*b1^3*mu^2+
          6*d(0,2)*b1*b2*mu+
          f(1,4)*b1^3*(1+3*mu+3*mu^2)+
          f(1,3)*b1*b2*(3+6*mu)+
          f(1,2)*b3+
          f(0,3)*b1^3*(3+6*mu)+
          6*f(0,2)*b1*b2
      )*a1*X[,vetor_theta[ordem_indices[1]]]
      )
    }
    if(nphi==4){ #k_phiphiphiphi
      kappa=sum(
        j(4)*b1^4+6*j(3)*b1^2*b2+j(2)*(3*b2^2+4*b1*b3)
      )
    }
    return(kappa)
  }
  
  k3 = function(r,s,t){
    vetor_theta = c(r,s,t) # vetor de índices do theta 0=phi
    ordem_indices = order(vetor_theta,decreasing=T)
    nphi=manyzeros(vetor_theta) # contando os phis
    
    if(nphi==0){ # k_rst
      kappa=sum((
        d(3,3)*a1^3+3*d(2,2)*a1*a2
      )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
      )
    }
    if(nphi==1){ #k_rsphi
      kappa=sum((
        d(2,3)*mu*a1^2+
          d(1,2)*(mu*a2+2*a1^2)+
          f(2,3)*a1^2+
          f(1,2)*a2
      )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
      )
    }
    if(nphi==2){ #k_rphiphi
      kappa=sum((
        d(1,3)*mu^2*b1^2+
          d(1,2)*mu*b2+
          2*d(0,2)*mu*b1^2+
          f(1,3)*(1+2*mu)*b1^2+
          f(1,2)*b2+
          2*f(0,2)*b1^2
      )*a1*X[,vetor_theta[ordem_indices[1]]]
      )
    }
    if(nphi==3){ #k_phiphiphi
      kappa=sum(
        j(3)*b1^3+3*j(2)*b1*b2
      )
    }
    return(kappa)
  }
  
  k31 = function(r,s,t,u){# k_rstu e afins
    vetor_theta = c(r,s,t,u) # vetor de índices do theta 0=phi
    ordem_indices = order(vetor_theta,decreasing=T)
    nphi=manyzeros(vetor_theta) # contando os phis
    
    if(u==0){
      if(nphi==1){# k_rst(phi)
        kappa=sum((
          d(3,4)*mu*a1^3+
            3*d(2,3)*a1*(a1^2+mu*a2)+
            6*d(1,2)*a1*a2+
            f(3,4)*a1^3+
            3*f(2,3)*a1*a2
        )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
        )
      } 
      if(nphi==2){# k_rsphi(phi)
        kappa=sum((
          d(2,4)*mu^2*a1^2*b1^2+
            d(2,3)*mu*a1^2*b2+
            d(1,3)*mu*b1^2*(4*a1^2+mu*a2)+
            d(1,2)*b2*(mu*a2+2*a1^2)+
            d(0,2)*b1^2*(mu*a2+2*a1^2)+
            f(2,4)*(1+2*mu)*a1^2*b1^2+
            f(2,3)*a1^2*b2+
            f(1,3)*b1^2*(4*a1^2+(1+2*mu)*a2)+
            f(1,2)*a2*b2+
            f(0,2)*a2*b1^2
        )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
        )
      }
      if(nphi==3){# k_rphiphi(phi)
        kappa=sum((
          d(1,4)*b1^3*mu^3+
            3*d(1,3)*b1*b2*mu^2+
            d(1,2)*b3*mu+
            3*d(0,3)*b1^3*mu^2+
            5*d(0,2)*b1*b2*mu+
            f(1,4)*b1^3*(1+3*mu+3*mu^2)+
            f(1,3)*b1*b2*(3+6*mu)+
            f(1,2)*b3+
            f(0,3)*b1^3*(3+6*mu)+
            5*f(0,2)*b1*b2
        )*a1*X[,vetor_theta[ordem_indices[1]]]
        )
      }
      if(nphi==4){# k_phiphiphi(phi)
        kappa=sum(
          j(4)*b1^4+6*j(3)*b1^2*b2+3*j(2)*(b2^2+b1*b3)
        )
      }
    }else{
      if(nphi==0){#k_rst(u)
        kappa=sum((
          d(4,4)*a1^4+6*d(3,3)*a1^2*a2+3*d(2,2)*(a2^2+a1*a3)
        )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]*X[,vetor_theta[ordem_indices[4]]]
        )
      }
      if(nphi==1){#k_rsphi(t)
        kappa=sum((
          d(3,4)*mu*a1^3+
            3*d(2,3)*(mu*a1*a2+a1^3)+
            d(1,2)*(mu*a3+5*a1*a2)+
            f(3,4)*a1^3+
            3*f(2,3)*a1*a2+
            f(1,2)*a3
        )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
        )
      }
      if(nphi==2){#k_rphiphi(s)
        kappa=sum((
          d(2,4)*mu^2*a1^2*b1^2+
            d(2,3)*mu*a1^2*b2+
            d(1,3)*mu*b1^2*(4*a1^2+mu*a2)+
            d(1,2)*b2*(mu*a2+a1^2)+
            2*d(0,2)*b1^2*(mu*a2+a1^2)+
            f(2,4)*(1+2*mu)*a1^2*b1^2+
            f(2,3)*a1^2*b2+
            f(1,3)*b1^2*(4*a1^2+(1+2*mu)*a2)+
            f(1,2)*a2*b2+
            2*f(0,2)*a2*b1^2
        )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
        )
      }
      if(nphi==3){# k_phiphiphi(r)
        kappa=sum((
          d(1,4)*mu^3*b1^3+
            3*d(1,3)*mu^2*b1*b2+
            3*d(0,3)*mu^2*b1^3+
            6*d(0,2)*mu*b1*b2+
            f(1,4)*(1+3*mu+3*mu^2)*b1^3+
            3*f(1,3)*(1+2*mu)*b1*b2+
            3*f(0,3)*(1+2*mu)*b1^3+
            6*f(0,2)*b1*b2
        )*a1*X[,vetor_theta[ordem_indices[1]]]
        )
      }
    }
    return(kappa)
  }
  
  k22 = function(r,s,t,u){
    vetor_theta = c(r,s,t,u) # vetor de índices do theta 0=phi
    ordem_indices = order(vetor_theta,decreasing=T)
    nphi=manyzeros(vetor_theta) # contando os phis
    
    if(u==0){
      if(t==0){
        if(nphi==2){# k_rs(phiphi)
          kappa=sum((
            d(2,4)*mu^2*b1^2+
              d(2,3)*mu*b2+
              4*d(1,3)*mu*b1^2+
              2*d(1,2)*b2+
              2*d(0,2)*b1^2+
              f(2,4)*(1+2*mu)*b1^2+
              f(2,3)*b2+
              4*f(1,3)*b1^2
          )*a1^2*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
          )
        }
        if(nphi==3){# k_rphi(phiphi)
          kappa=sum((
            d(1,4)*b1^3*mu^3+
              3*d(1,3)*b1*b2*mu^2+
              d(1,2)*b3*mu+
              2*d(0,3)*b1^3*mu^2+
              3*d(0,2)*b1*b2*mu+
              f(1,4)*b1^3*(1+3*mu+3*mu^2)+
              f(1,3)*b1*b2*(3+6*mu)+
              f(1,2)*b3+
              f(0,3)*b1^3*(2+4*mu)+
              3*f(0,2)*b1*b2
          )*a1*X[,vetor_theta[ordem_indices[1]]]
          )
        }
        if(nphi==4){#k_phiphi(phiphi)
          kappa=sum(
            j(4)*b1^4+5*j(3)*b1^2*b2+2*j(2)*(b2^2+b1*b3)
          )
        }
      }else{
        if(nphi==1){#k_rs(tphi)
          kappa=sum((
            d(3,4)*mu*a1^3+
              d(2,3)*a1*(3*a1^2+2*mu*a2)+
              4*d(1,2)*a1*a2+
              f(3,4)*a1^3+
              2*f(2,3)*a1*a2
          )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
          )
        }
        if(nphi==2){#k_rphi(sphi)
          kappa=sum((
            d(2,4)*mu^2*a1^2*b1^2+
              d(2,3)*mu*a1^2*b2+
              d(1,3)*mu*b1^2*(3*a1^2+mu*a2)+
              d(1,2)*b2*(mu*a2+a1^2)+
              d(0,2)*b1^2*(mu*a2+a1^2)+
              f(2,4)*(1+2*mu)*a1^2*b1^2+
              f(2,3)*a1^2*b2+
              f(1,3)*b1^2*(3*a1^2+(1+2*mu)*a2)+
              f(1,2)*a2*b2+
              f(0,2)*a2*b1^2
          )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
          )
        }
        if(nphi==3){#k_phiphi(rphi)
          kappa=sum((
            d(1,4)*mu^3*b1^3+
              2*d(1,3)*mu^2*b1*b2+
              3*d(0,3)*mu^2*b1^3+
              4*d(0,2)*mu*b1*b2+
              f(1,4)*(1+3*mu+3*mu^2)*b1^3+
              2*f(1,3)*(1+2*mu)*b1*b2+
              3*f(0,3)*(1+2*mu)*b1^3+
              4*f(0,2)*b1*b2
          )*a1*X[,vetor_theta[ordem_indices[1]]]
          )
        }
      }
    }else{
      if(t==0){
        if(nphi==1){#k_rs(phit)
          kappa=sum((
            d(3,4)*mu*a1^3+
              d(2,3)*a1*(3*a1^2+2*mu*a2)+
              4*d(1,2)*a1*a2+
              f(3,4)*a1^3+
              2*f(2,3)*a1*a2
          )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
          )
        }
        if(nphi==2){#k_rphi(phis)
          kappa=sum((
            d(2,4)*mu^2*a1^2*b1^2+
              d(2,3)*mu*a1^2*b2+
              d(1,3)*mu*b1^2*(3*a1^2+mu*a2)+
              d(1,2)*b2*(mu*a2+a1^2)+
              d(0,2)*b1^2*(mu*a2+a1^2)+
              f(2,4)*(1+2*mu)*a1^2*b1^2+
              f(2,3)*a1^2*b2+
              f(1,3)*b1^2*(3*a1^2+(1+2*mu)*a2)+
              f(1,2)*a2*b2+
              f(0,2)*a2*b1^2
          )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
          )
        }
        if(nphi==3){#k_phiphi(phir)
          kappa=sum((
            d(1,4)*mu^3*b1^3+
              2*d(1,3)*mu^2*b1*b2+
              3*d(0,3)*mu^2*b1^3+
              4*d(0,2)*mu*b1*b2+
              f(1,4)*(1+3*mu+3*mu^2)*b1^3+
              2*f(1,3)*(1+2*mu)*b1*b2+
              3*f(0,3)*(1+2*mu)*b1^3+
              4*f(0,2)*b1*b2
          )*a1*X[,vetor_theta[ordem_indices[1]]]
          )
        }
      }else{
        if(nphi==0){#k_rs(tu)
          kappa=sum((
            d(4,4)*a1^4+5*d(3,3)*a1^2*a2+2*d(2,2)*(a2^2+a1*a3)
          )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]*X[,vetor_theta[ordem_indices[4]]]
          )
        }
        if(nphi==1){#k_rphi(st)
          kappa=sum((
            d(3,4)*mu*a1^3+
              d(2,3)*(3*mu*a1*a2+2*a1^3)+
              d(1,2)*(mu*a3+3*a1*a2)+
              f(3,4)*a1^3+
              3*f(2,3)*a1*a2+
              f(1,2)*a3
          )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
          )
        }
        if(nphi==2){#k_phiphi(rs)
          kappa=sum((
            d(2,4)*mu^2*a1^2+
              d(1,3)*mu*(4*a1^2+mu*a2)+
              2*d(0,2)*(mu*a2+a1^2)+
              f(2,4)*(1+2*mu)*a1^2+
              f(1,3)*(4*a1^2+(1+2*mu)*a2)+
              2*f(0,2)*a2
          )*b1^2*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
          )
        }
      }
    }
    return(kappa)
  }
  
  k21 = function(r,s,t){
    vetor_theta = c(r,s,t) # vetor de índices do theta 0=phi
    ordem_indices = order(vetor_theta,decreasing=T)
    nphi=manyzeros(vetor_theta) # contando os phis
    if(t==0){
      if(nphi==1){#k_rs(phi)
        kappa=sum((
          d(2,3)*mu+2*d(1,2)+f(2,3)
        )*b1*a1^2*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]])
      }
      if(nphi==2){#k_rphi(phi)
        kappa=sum((
          d(1,3)*mu^2*b1^2+
            d(1,2)*mu*b2+
            d(0,2)*mu*b1^2+
            f(1,3)*(1+2*mu)*b1^2+
            f(1,2)*b2+
            f(0,2)*b1^2
        )*a1*X[,vetor_theta[ordem_indices[1]]]
        )
      }
      if(nphi==3){#k_phiphi(phi)
        kappa=sum(
          j(3)*b1^3+2*j(2)*b1*b2
        )
      }
    }else{
      if(nphi==0){#k_rs(t)
        kappa=sum((
          d(3,3)*a1^3+2*d(2,2)*a1*a2
        )*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]*X[,vetor_theta[ordem_indices[3]]]
        )
      }
      if(nphi==1){#k_rphi(s)
        kappa=sum((
          d(2,3)*mu*a1^2+
            d(1,2)*(mu*a2+a1^2)+
            f(2,3)*a1^2+
            f(1,2)*a2
        )*b1*X[,vetor_theta[ordem_indices[1]]]*X[,vetor_theta[ordem_indices[2]]]
        )
      }
      if(nphi==2){#k_phiphi(r)
        kappa=sum((
          d(1,3)*mu^2+
            2*d(0,2)*mu+
            f(1,3)*(1+2*mu)+
            2*f(0,2)
        )*a1*b1^2*X[,vetor_theta[ordem_indices[1]]])
      }
    }
    return(kappa)
  }
  
  
  
  #### matrizes ####
  A=array(rep(0,k^4),dim=c(k,k,k,k))
  P=array(rep(0,k^3),dim=c(k,k,k))
  Q=P
  for(r in 1:k){
    ri=ifelse(r>p,0,r)
    for(s in 1:k){
      si=ifelse(s>p,0,s)
      for(t in 1:k){
        ti=ifelse(t>p,0,t)
        P[t,r,s]=k3(ri,si,ti)
        for(u in 1:k){
          ui=ifelse(u>p,0,u)
          A[t,u,r,s]=k4(ri,si,ti,ui)/4-k31(ri,si,ti,ui)+k22(ri,ti,si,ui)
        }
      }
    }
  }
  
  for(u in 1:k){
    ui=ifelse(u>p,0,u)
    for(r in 1:k){
      ri=ifelse(r>p,0,r)
      for(s in 1:k){
        si=ifelse(s>p,0,s)
        Q[u,r,s]<-k21(si,ui,ri)
      }
    }
  }
  
  
  L=matrix(rep(0,k^2),ncol=k)
  M1<-M2<-M3<-N1<-N2<-N3<-L
  for(r in 1:k){
    for(s in 1:k){
      L[r,s]=sum(diag(K_1%*%A[r,s,,]))
      M1[r,s]=sum(diag(K_1%*%P[r,,]%*%K_1%*%P[s,,]))
      M2[r,s]=sum(diag(K_1%*%P[r,,]%*%K_1%*%t(Q[s,,])))
      M3[r,s]=sum(diag(K_1%*%Q[r,,]%*%K_1%*%Q[s,,]))
      N1[r,s]=sum(diag(P[r,,]%*%K_1))*sum(diag(P[s,,]%*%K_1))
      N2[r,s]=sum(diag(P[r,,]%*%K_1))*sum(diag(Q[s,,]%*%K_1))
      N3[r,s]=sum(diag(Q[r,,]%*%K_1))*sum(diag(Q[s,,]%*%K_1))
    }
  }
  M=-M1/6+M2-M3
  N=-N1/4+N2-N3
  return(sum(diag(K_1%*%(L-M-N))))
}
