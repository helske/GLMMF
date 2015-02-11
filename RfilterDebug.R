
r<-matrix(0,14,12)
t<-28
for(i in 12:2)
  r[,i-1]<-o0$model$Z[i,,t]*o0$v[t,i]/o0$F[i,t]+ 
  diag(14)%*%r[,i]-((o0$model$Z[i,,t])%*%t(o0$K[,i,t])/o0$F[i,t])%*%r[,i]

t<
all.equal(((o0$model$Z[i,,t])%*%t(o0$K[,i,t])/o0$F[i,t])%*%r[,i],
          ((matrix(o0$model$Z[i,c(2,13,14),t],ncol=1)%*%t(o0$K[,i,t])/o0$F[i,t])%*%r[,i])
          )


(o0$K[,i,t]%*%t(o0$model$Z[i,c(2,13,14),t])/o0$F[i,t])%*%r[c(2,13,14),i]

((o0$K[,i,t])%*%t(o0$model$Z[i,,t])/o0$F[i,t])[,c(2,13,14)]
##

load("kmodel.rda")
library(KFAS)
app<-approxSSM(kmodel)
out<-KFS(app,simpl=F,filtering="state")

vt<-Ft<-matrix(0,28,12)
Kt<-array(0,c(14,12,28))
Pt<-array(0,c(14,14,13,29))
at<-array(0,c(14,13,29))
yt<-app$y
Z<-app$Z
Tt<-app$T[,,1]
Qt<-diag(c(rep(0,12),1,1))
Pt[,,1,1]<-diag(c(rep(1e7,12),1,1))

for(t in 1:28){
  for(i in 1:12){
    vt[t,i]<-yt[t,i]-Z[i,,t]%*%at[,i,t]
    Ft[t,i]<-t(Z[i,,t])%*%Pt[,,i,t]%*%(Z[i,,t]) + app$H[i,i,t]
    Kt[,i,t]<- Pt[,,i,t]%*%(Z[i,,t])
    
    at[,i+1,t]<-at[,i,t]+Kt[,i,t]/Ft[t,i]*vt[t,i]
    Pt[,,i+1,t]<-Pt[,,i,t] - Kt[,i,t]%*%t(Kt[,i,t])/Ft[t,i]   
    
    
  }

    at[,1,t+1]<-Tt%*%(at[,13,t])
    Pt[,,1,t+1]<-Tt%*%(Pt[,,13,t])%*%t(Tt) + Qt 

}

all.equal(out$a[29,1:12],out$alpha[28,1:12])
all.equal(out$a[29,],at[,1,29],check.attributes=FALSE)
all.equal(out$a[20,],at[,1,20],check.attributes=FALSE)

app2<-app
app2$a1[1:12]<-out$alpha[28,1:12]
app2$P1inf[]<-0
out2<-KFS(app2,simpl=F,filtering="state")

all.equal(out$alpha,out2$alpha) #true
out2$K
all.equal(t(app$Z[,,28])%*%solve(app$Z[,,28]%*%out2$P[,,28]%*%t(app$Z[,,28])+app$H[,,28])%*%app$Z[,,28],out2$N[,,28],
          check.attributes=FALSE)

t<-10
Tt<-app$T[,,1]
Finv<-solve(app$Z[,,t]%*%out2$P[,,t]%*%t(app$Z[,,t])+app$H[,,t])
L<-diag(14)-out2$K[,,t]%*%Finv%*%app$Z[,,t] #toinen termi on nolla ei-faktoreille
t(Tt)%*%out2$N[,,t]%*%Tt #nolla faktoreille
all.equal(t(L)%*%t(Tt)%*%out2$N[,,t]%*%Tt%*%L,t(Tt)%*%out2$N[,,t]%*%Tt)

all.equal(t(L)%*%t(Tt)%*%out2$N[,,t]%*%Tt%*%L)
all.equal(out2$N[,,t-1],t(app$Z[,,t])%*%solve(app$Z[,,t]%*%out2$P[,,t]%*%t(app$Z[,,t])+app$H[,,t])%*%app$Z[,,t]+
            t(L)%*%t(Tt)%*%out2$N[,,t]%*%Tt%*%L,t(Tt)%*%out2$N[,,t]%*%Tt)

t<-20
Finv<-solve(app$Z[,,t]%*%out2$P[,,t]%*%t(app$Z[,,t])+app$H[,,t])
L<-diag(14)-out2$K[,,t]%*%Finv%*%app$Z[,,t]
all.equal(t(app$Z[,,t])%*%solve(app$Z[,,t]%*%out2$P[,,t]%*%t(app$Z[,,t])+app$H[,,t])%*%app$Z[,,t]+
            t(L)%*%t(Tt)%*%out2$N[,,t+1]%*%Tt%*%L,out2$N[,,t],
          check.attributes=FALSE)

Finv<-solve(app$Z[,,t]%*%out2$P[,,t]%*%t(app$Z[,,t])+app$H[,,t])
all.equal(t(app$Z[,,t])%*%Finv%*%app$Z[,,t]+
            t(Tt)%*%out2$N[,,t+1]%*%Tt,out2$N[,,t],
          check.attributes=FALSE) ## EI TARVITA K:TA, PELKKÄ F ja P riittää!!
all.equal(out2$N[,,t],out$N[,,t]) #NOT TRUE, eli tarvitaan kuitenkin:
Finv<-solve(app$Z[,,t]%*%out$P[,,t]%*%t(app$Z[,,t])+app$H[,,t])
all.equal(t(app$Z[,,t])%*%Finv%*%app$Z[,,t]+
            t(Tt)%*%out$N[,,t+1]%*%Tt,out$N[,,t],
          check.attributes=FALSE)

Finv<-solve(app$Z[,,t]%*%out$P[,,t]%*%t(app$Z[,,t])+app$H[,,t])
L<-diag(14)-out$K[,,t]%*%Finv%*%app$Z[,,t]
r<-t(app$Z[,,t])%*%Finv%*%out$v[t,]+t(L)%*%t(Tt)%*%out$r[,t+1]
out$r[,t]
c(r)
out$alpha[t+1,]
Tt%*%out$alpha[t,]+app$R[,,1]%*%diag(2)%*%t(app$R[,,1])%*%out$r[,t+1]

###
L<-diag(14)-out$K[,,t]%*%Finv%*%app$Z[,,t]
r<-t(app$Z[,,t])%*%Finv%*%out$v[t,]+t(L)%*%t(Tt)%*%out$r[,t+1]
out$r[,t]
c(r)
out$alpha[t+1,]
Tt%*%out$alpha[t,]+app$R[,,1]%*%diag(2)%*%t(app$R[,,1])%*%out$r[,t+1]
##
all.equal(out$F,out2$F)
