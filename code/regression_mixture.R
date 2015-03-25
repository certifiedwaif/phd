# regression_mixture.R

# @return A list of generated data
generate <- function()
{
	n <- 30
	mX <- cbind(rep(1, n), runif(n, -1, 1))
	alpha_true_model <- c(1, -1)
	beta_true_model <- c(-1, 1)
	prob <- 0.5
	A_1 <- 0.1
	A_2 <- 0.1
	B <- 0.1
	sigma2_alpha <- 1.0
	sigma2_beta <- 1.5
	vy <- rep(NA, n)
	class_labels <- rep(NA, n)
	
	for (i in 1:n) {
		class_labels[i] <- rbinom(1, 1, prob)
		if (class_labels[i] <- 0) {
			vy[i] <- mX[i, ] %*% alpha_true_model + rnorm(1, 0, sigma2_alpha)
		} else {
			vy[i] <- mX[i, ] %*% beta_true_model + rnorm(1, 0, sigma2_beta)
		}
	}

	return(list(A_1=A_1, A_2=A_2, B=B,
				n=n, mX=mX, alpha_true_model=alpha_true_model,
				beta_true_model=beta_true_model, prob=prob, vy=vy))
}

# @param data - The data list returned by generate
# @return A list of the fitted parameters
vb <- function(data)
{
	n = data$n
	mX = data$mX
	vy = data$vy
	A_1 = data$A_1
	A_2 = data$A_2
	B = data$B

	# Initialise
	ITERATIONS <- 10
	p <- 2
	vzero <- rep(0, n)
	vhalf <- rep(0.5, n)
	mI <- diag(1, p)
	m_S_alpha <- mI
	m_S_beta <- mI
	vm_alpha <- vzero
	vm_beta <- vzero
	veta <- vzero
	vp <- vhalf
	mP <- diag(vp)
	r_1 <- A_1
	psi_1 <- B
	r_0 <- A_2
	psi_0 <- B
	sigmainv2_alpha <- r_1/psi_1
	sigmainv2_beta <- r_0/psi_0
  
	for (i in 1:ITERATIONS) {
		# Mean field update
		m_S_alpha <- solve(r_1/psi_1)*t(mX)%*%mP%*%mX + sigmainv2_alpha * mI
		vm_alpha <- m_S_alpha %*% t(mX) %*% mP %*% vy * r_1/psi_1
		m_S_beta <- solve(r_0/psi_0)*t(mX)%*%mP%*%mX + sigmainv2_beta * mI
		vm_beta <- m_S_beta %*% t(mX) %*% mP %*% vy * r_0/psi_0
		veta <- logit(rho) + (-0.5*log(2*pi) - 0.5*log(psi_1) + 0.5*digamma(r_1) - 0.5*r_1/psi_1*((vy - t(mX)%*%vm_alpha)^2 + t(mX)%*%mS_alpha%*%mX))
		veta <- veta - (-0.5*log(2*pi) - 0.5*log(psi_0) + 0.5*digamma(r_0) - 0.5*r_0/psi_0*((vy - t(mX)%*%vm_beta)^2 + t(mX)%*%mS_beta%*%mX))
		vp <- expit(veta)
		mP <- diag(vp)
		r_1 <- A_1 + 0.5 * sum(vp)
		psi_1 <- B + 0.5 * (t(vp) %*% (vy - mX%*%vm_alpha)^2 + tr(m_S_alpha %*% t(mX) %*% mP %*% mX))
		r_0 <- A_2 + 0.5 * sum(vp)
		psi_0 <- B + 0.5 * (t(vp) %*% (vy - mX%*%vm_beta)^2 + tr(m_S_beta %*% t(mX) %*% mP %*% mX))
	}
	return(list(m_S_alpha <- m_S_alpha, m_S_beta <- m_S_beta, vm_alpha <- vm_alpha,
				vm_beta <- vm_beta, veta <- veta, vp <- vp, mP <- mP, r_1 <- r_1,
				psi_1 <- psi_1, r_0 <- r_0, psi_0 <- psi_0))
}

main <- function()
{
	data <- generate()
	result <- vb(data)
	print(result)
}
main()