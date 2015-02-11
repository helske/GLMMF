## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
d.AD<-rbind(d.AD,d.AD,d.AD,d.AD,d.AD,d.AD,d.AD,d.AD,d.AD)
d.AD<-rbind(d.AD,d.AD,d.AD,d.AD,d.AD,d.AD,d.AD,d.AD,d.AD)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
anova(glm.D93)
summary(glm.D93)

model<-buildSSGLMM(counts~outcome+treatment,data=d.AD,dist="poisson")
dist<-pmatch(x = model$distribution, 
             table = c("gaussian", "poisson", "binomial", 
                       "gamma", "negative binomial"), duplicates.ok = TRUE)
theta<-KFAS:::init_theta(model$y, model$u, model$distribution)

f1<-function() {
  

  
 
  c(expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                     dist, model$tol,100,1e-8,theta)$theta)
}
f2<-function() c(glm(counts ~ outcome + treatment, data=d.AD,family = poisson(),etas=theta)$line)
library(microbenchmark)
all.equal(f1(),f2(),check.attributes=FALSE)
microbenchmark(f1(),f2())

## an example with offsets from Venables & Ripley (2002, p.189)
utils::data(anorexia, package = "MASS")

anorex.1 <- glm(Postwt ~ Prewt + Treat + (Prewt),
                family = gaussian, data = anorexia)
summary(anorex.1)

########

model<-buildSSGLMM(counts~outcome+treatment,data=d.AD,dist="poisson")
dist<-pmatch(x = model$distribution, 
             table = c("gaussian", "poisson", "binomial", 
                       "gamma", "negative binomial"), duplicates.ok = TRUE)
theta<-KFAS:::init_theta(model$y, model$u, model$distribution)

f1<-function() {
  
  
  
  
  c(expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                    dist, model$tol,100,1e-8,theta)$theta)
}
f2<-function() c(glm(counts ~ outcome + treatment, data=d.AD,family = poisson(),etas=theta)$line)
library(microbenchmark)
all.equal(f1(),f2(),check.attributes=FALSE)
microbenchmark(f1(),f2())
