##
## ISO 28037:2010
##

###############################
## Section 10. x,Ux,y,Uy,Uxy ##
###############################

ggmr <- function(x, y, Ux = diag(0, length(x)), 
			Uy = diag(1, length(x)), 
			Uxy = diag(0, length(x)),
			subset = rep(TRUE, length(x)), 
			tol = sqrt(.Machine$double.eps), 
			max.iter = 100, 
			alpha = 0.05,
			coef.H0 = c(0, 1)) {
	x <- x[subset]
	y <- y[subset]
	Ux <- Ux[subset,subset]
	Uy <- Uy[subset,subset]
	Uxy <- Uxy[subset,subset]

	m <- length(x)
	## initial validation
	if (length(y) != m) 
		stop("response and predictor vectors must have same length.")
	if (any(dim(Ux) != c(m, m))) stop(paste("the dimension of the predictor ",
		"covariance matrix must equate the predictor vector length."))
	if (any(dim(Uy) != c(m, m))) stop(paste("the dimension of the response ",
		"covariance matrix must equate the response vector length."))
	if (any(dim(Uxy) != c(m, m))) stop(paste("the dimension of the ",
		"predictor-response covariance matrix must equate the predictor ",
		"vector length."))

#	if ((sum(diag(Ux)==0)==m) || (sum(diag(Ux)!=0)!=m)) stop()
#	if (det(Ux)==0)

	if (max.iter > 1000) max.iter <- 1000
	if (max.iter < 10) max.iter <- 10
	
	if (all(Ux==0) && all(Uxy==0) && all(diag(Uy)>0) && all((Uy-diag(diag(Uy)))==0)) {
		## ISO 28037, section 6
		X<- matrix(cbind(rep(1, m), x), m, 2)
		W<- Uy
		diag(W)<- 1/diag(Uy)
		Sigma.betas<- ginv(t(X)%*%W%*%X)
		betas<- Sigma.betas%*%t(X)%*%W%*%y
		xi <- rep(NA, m)
		y.hat <- X%*%betas
		residuals <- y-y.hat
		test.stat.validation <- sum(diag(W)*residuals^2)
		test.stat <- t(c(betas-coef.H0))%*%ginv(Sigma.betas)%*%c(betas-coef.H0)
		p.value <- 1 - pchisq(test.stat.validation, m - 2)
		ii <- 0
		curr.tol <- 0
	} else if (all(diag(Ux)>0) && all(t(Ux)==Ux) && all(diag(Uy)>0) && all(t(Uy)==Uy) && all(t(Uxy)==Uxy)) {
		if (any(eigen(Ux)$values<=0)) stop("Ux is expected to be a positive definite symmetric matrix.")
		if (any(eigen(Uy)$values<=0)) stop("Uy is expected to be a positive definite symmetric matrix.")
		## ISO 28037, section 7 & 10
		U <- diag(0, 2*m, 2*m)
		U[1:m, 1:m] <- Ux
		U[(m + 1):(2*m), (m + 1):(2*m)] <- Uy
		U[(m + 1):(2*m), 1:m] <- Uxy
		U[1:m, (m + 1):(2*m)] <- t(Uxy)
		if (any(eigen(U)$values<=0)) stop("Combined covariance matrix is expected to be a positive definite symmetric matrix.")
		if (det(U) == 0) stop("Combined covariance matrix is non-invertible.")
		## step 3
		L <- t(chol(U))
		## step 1
		t <- c(x, lm(y~x, weights = 1/diag(Uy))$coef[1:2])
		names(t) <- NULL
		dt <- rep(Inf, m + 2)
		ii <- 1
		tt <- matrix(NA, max.iter, m + 2)
		tt[ii, ] <- t
		## step 10
		f.tilde <- rep(NA, length(x)+length(y))
		while (any(dt > tol*t) & (ii < max.iter)) {
			## step 2
			f <- c(x, y) - c(t[1:m], t[m + 1] + t[m + 2]*t[1:m])
			J <- matrix(0, 2*m, m + 2)
			J[1:m, 1:m] <- -diag(1, m)
			J[(m + 1):(2*m), 1:m] <- -t[m + 2]*diag(1, m)
			J[(m + 1):(2*m), m + 1] <- -1
			J[(m + 1):(2*m), m + 2] <- -t[1:m]
			## step 4
			f.tilde <- ginv(t(L) %*% L) %*% t(L) %*% f
			J.tilde <- ginv(t(L) %*% L) %*% t(L) %*% J
			## step 5
			g <- t(J.tilde) %*% f.tilde
			H <- t(J.tilde) %*% J.tilde
			## step 6
			M <- t(chol(H))
			## step 7
			q <- -ginv(t(M) %*% M) %*% t(M) %*% g
			## step 8
			dt <- ginv(M %*% t(M)) %*% M %*% q
			## step 9
			t <- t + dt
			ii <- ii + 1
			tt[ii, ] <- t
		}
		## step 11
		xi <- t[1:m]
		betas <- t[c(m + 1,m + 2)]
		mm <- M[c(m + 1,m + 2), c(m + 1,m + 2)]
		Sigma.betas <- mm
		Sigma.betas[1, 1] <- (mm[2,1]^2 + mm[2,2]^2)/(prod(diag(mm))^2)
		Sigma.betas[2, 2] <- 1/(mm[2,2]^2)
		Sigma.betas[1, 2] <- Sigma.betas[2,1] <- -mm[2,1]/(mm[1,1]*mm[2,2]^2)

		## hypothesis testing againt a specified coef.H0
		delta <- (betas - coef.H0)
		test.stat <- delta %*% ginv(Sigma.betas) %*% delta

		p.value <- 1 - pchisq(sum(f.tilde^2), m - 2)
		curr.tol <- max(dt/t)

		## model validation test
		if (all((Ux - diag(Ux))== 0) && all((Uy-diag(Uy))==0) && all(Uxy==0)) {
			# for section 7
			ti <- 1/(diag(Uy)+betas[2]^2*diag(Ux))
			zi <- y-betas[1]-betas[2]*x
			fi <- sqrt(ti)
			gi <- fi*xi
			hi <- fi*zi
			g0 <- sum(fi*gi)/sum(fi^2)
			h0 <- sum(fi*hi)/sum(fi^2)
			g.tilde <- gi - g0*fi
			h.tilde <- hi - h0*fi
			G2.tilde <- sum(g.tilde^2)
			delta.b <- sum(g.tilde*h.tilde)/G2.tilde
			ri <- h.tilde - delta.b*g.tilde
			test.stat.validation <- sum(ri^2)
		} else {
			## for section 10.
			test.stat.validation <- sum(f.tilde^2)
		}
	} else {
		stop(paste("symetric matrices are expected,",
			"mixed constant and random variables for the predictor or response are not allowed."))
	}

	res <- list(coefficients = betas, 
			cov = Sigma.betas, 
			xi = xi,
			chisq.validation = test.stat.validation,
			chisq.ht = test.stat,
			chisq.cri = qchisq(1 - alpha, m - 2),
			p.value = p.value,
			curr.iter = ii, curr.tol = curr.tol
		)
	class(res) <- "ggmr"
	return(res)
}
