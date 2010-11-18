inudge.fit <-
function(data, K = 2, weights = NULL, pi = NULL, mu = NULL,
sigma = NULL, tol=1e-16, max.iter=2000, z = NULL)
{
  
  x <- unlist(data);
  n <- length(x);
  if(is.null(weights)) weights <- rep(1,length(x));
         
	if(K < 1) stop("It is expected that there is at least one normal component; therefore, K needs to be greater than 0.");
  # Uniform-Normal^K model 
    a <- min(x);
    b <- max(x);
    r <- matrix(0,n,K);
    m <- mean(x);
    s <- sqrt(var(x));
    d <- abs((x-m)/s) > 2;
   if(is.null(z)){
    z <- matrix(0,n,2);
    z[,2] <- d^2;
    z[,1] <- 1 - z[,2];
   }
    
  # initialize probability of uniform component
  if(is.null(pi))  pi <- rep(sum(z[,1])/(n*K),K); 
  # initialize normal parameters
  if(is.null(mu)){ muhat <- runif(K,sum(z[, 1] * x)/sum(z[, 1])-s,sum(z[, 1] * x)/sum(z[, 1])+s);
  }else{ muhat <- mu;}
 
  if(is.null(sigma)){
    sigma2hat <- sample(sum(z[, 1] * (rep(x,K) - muhat)^2)/sum(z[, 1]),K,replace = TRUE);
    sigmahat <- sqrt(sigma2hat);
  }else{ sigmahat <- sigma;}
  
  llike <- c(0,100);
  criterion <- abs(llike[1]-llike[2]);
  iter <- 0;
  f1 <- matrix(0,n,1);
  f0 <- matrix(0,n,K);
  
  while((criterion > tol) & (iter < max.iter)){
     # save current estimation in case of non-convergence
     tmpmuhat <- muhat;
     tmpsigmahat <- sigmahat;
     tmppi <- pi;
                   
  	iter <- iter + 1;
  	# E-step
  	for (k in 1:K){
  		f0[,k] <- pi[k] * dnorm(x, muhat[k],sigmahat[k]);
  	}
  	f1 <- (1- sum(pi)) * dunif(x,a,b);
  	for (k in 1:K){
  	r[,k] <- pi[k] * dnorm(x, muhat[k],sigmahat[k]) /
  			 (apply(f0,1,sum)+ f1);
  	}
  	z[,2] <- f1 / (apply(f0,1,sum)+ f1);

  	# M-step
  	pi <- apply(r,2,sum)/n;
  	for (k in 1:K){
  		muhat[k] <- sum(r[,k] * x)/sum(r[,k]);
  		sigmahat[k] <- sqrt(sum(r[,k] * (x - muhat[k])^2)/sum(r[,k]));
  	}
  	# avoiding error due to no convergence (i.e. llike is NA)
  	tmpllike <- sum(log(apply(f0,1,sum)+f1));
  	if (!is.nan(tmpllike)&&!is.infinite(tmpllike)){
    	llike[2] <- llike[1];
    	llike[1] <- sum(log(apply(f0,1,sum)+f1));
    }else{
    muhat <- tmpmuhat;
    sigmahat <- tmpsigmahat;
    pi <- tmppi;
    }
  	criterion <- abs(llike[1] - llike[2]);
  }
  AIC <- llike[1] - (3*K);
  BIC <- 2*llike[1] - (3*K)*log(n);
  if (K ==1){
  model.name = "NUDGE";
  }else{
  model.name = "iNUDGE"
  }
  param_estim <- list(pdiff = z[, 2],mu=muhat,sigma=sigmahat,loglike=llike[1],AIC=AIC,BIC=BIC,
		iter=iter,range=c(a,b),a=a,b=b,K = K, pi = c(1-sum(pi),pi),phi=cbind(f0,f1),
    name=model.name);
	return(param_estim);
}

