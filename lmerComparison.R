library(lme4)
library(GLMMF)

(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
summary(fm1)# (with its own print method)

model<-  buildGLMMF(group="Subject",response="Reaction",common=~Days,random=~Days,data=sleepstudy,
                    random.cov=VarCorr(fm1)$Subject,u=attr(VarCorr(fm1), "sc")^2)

theta <- KFAS:::init_theta(model$y, model$u, rep(model$distribution,18 ))
dist<-rep(pmatch(x = model$distribution, 
             table = c("gaussian", "poisson", "binomial", 
                       "gamma", "negative binomial"), duplicates.ok = TRUE),18)

res<-GLMMF:::expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                             dist, model$tol, 100, 1e-10, theta, model$Zind, model$nfactors)


all.equal(res$logLik+KFAS:::scaling(model$y,model$u,rep(model$distribution,18),res$theta),
          as.numeric(logLik(fm1)),check.attributes=FALSE)

fit<-glmmf(group="Subject",response="Reaction",common.f=~Days,random=~Days,data=sleepstudy,
           ,control=list(trace=1,REPORT=1),init.disp=attr(VarCorr(fm1), "sc")^2,
           init.random.cov=VarCorr(fm1)$Subject,correlating.effects=TRUE)

fit<-glmmf(group="Subject",response="Reaction",common.f=~Days,random=~Days,data=sleepstudy,
           control=list(trace=1,REPORT=1),correlating.effects=TRUE)

fit<-glmmf(group="Subject",response="Reaction",common.f=~Days,random=~Days,data=sleepstudy,
           method="BFGS",control=list(trace=1,REPORT=1),correlating.effects=TRUE)

fit$model$u[1]
(attr(VarCorr(fm1), "sc")^2)

(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy,REML=TRUE))
logLik(fm1,REML=TRUE)
logLik(fm1,REML=FALSE)

(fm1 <- glmer(Reaction ~ Days + (Days | Subject), sleepstudy,REML=TRUE))
logLik(fm1,REML=TRUE)
logLik(fm1,REML=FALSE)
