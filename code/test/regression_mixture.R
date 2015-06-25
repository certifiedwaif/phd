# regression_mixture.R

# @return A list of generated data
generate <- function()
{
	n <- 100
	mX <- cbind(rep(1, n), runif(n, -1, 1))
	alpha_true_model <- c(1, -1)
	beta_true_model <- c(-1, 1)
	prob <- 0.5
	A_1 <- 1e5
	A_2 <- 1e5
	B <- 1e5
	sigma2_alpha <- 1.0
	sigma2_beta <- 1.0
	vy <- rep(NA, n)
	class_labels <- rbinom(n, 1, prob)

  for (i in 1:n) {
		if (class_labels[i] == 0) {
			vy[i] <- mX[i, ] %*% alpha_true_model + rnorm(1, 0, sqrt(sigma2_alpha))
		} else {
			vy[i] <- mX[i, ] %*% beta_true_model + rnorm(1, 0, sqrt(sigma2_beta))
		}
	}

	return(list(A_1=A_1, A_2=A_2, B=B,
							n=n, mX=mX, alpha_true_model=alpha_true_model,
							beta_true_model=beta_true_model, prob=prob, vy=vy))
}

logit <- function(p)
{
	log(p/(1-p))
}

expit <- function(x)
{
	exp(x)/(1 + exp(x))
}

tr <- function(m) sum(diag(m))

# @param data - The data list returned by generate
# @return A list of the fitted parameters
vb <- function(data)
{
	n <- data$n
	mX <- data$mX
	vy <- data$vy
	A_1 <- data$A_1
	A_2 <- data$A_2
	B <- data$B

	# Initialise
	ITERATIONS <- 10
	p <- 2
	vzero <- rep(0, n)
	vhalf <- rep(0.5, n)
	mI_p <- diag(1, p)
  mI_n <- diag(1, n)
	mS_alpha <- mI_p
	mS_beta <- mI_p
	vm_alpha <- c(1, 1)
	vm_beta <- c(1, -1)
	veta <- vzero
	vp <- rbinom(n, 1, 0.5)
	mP <- diag(vp, n, n)
	rho <- 0.5
	r_1 <- A_1
	psi_1 <- B
	r_0 <- A_2
	psi_0 <- B
	sigmainv2_alpha <- r_1/psi_1
	sigmainv2_beta <- r_0/psi_0
  
	for (i in 1:ITERATIONS) {
		# Mean field update
		mS_alpha <- solve((r_1/psi_1) * t(mX) %*% mP %*% mX + sigmainv2_alpha * mI_p)
		vm_alpha <- mS_alpha %*% t(mX) %*% mP %*% vy * r_1/psi_1
		mS_beta <- solve((r_0/psi_0) * t(mX) %*% (mI_n - mP) %*% mX + sigmainv2_beta * mI_p)
		vm_beta <- mS_beta %*% t(mX) %*% (mI_n - mP) %*% vy * r_0/psi_0
		veta <- logit(rho) + (-0.5 * log(2 * pi) - 0.5 * log(psi_1) + 0.5 * digamma(r_1) - 0.5 * r_1/psi_1*((vy - mX %*% vm_alpha)^2 + tr(mX %*% mS_alpha %*% t(mX))))
		veta <- veta - (-0.5 * log(2 * pi) - 0.5 * log(psi_0) + 0.5 * digamma(r_0) - 0.5 * r_0/psi_0*((vy - mX %*% vm_beta)^2 + tr(mX %*% mS_beta %*% t(mX))))
		vp <- expit(veta)
		diag(mP) <- vp
		r_1 <- A_1 + 0.5 * sum(vp)
		psi_1 <- as.vector(B + 0.5 * (t(vp) %*% (vy - mX %*% vm_alpha)^2 + tr(mS_alpha %*% (t(mX) %*% mP %*% mX))))
		r_0 <- A_2 + 0.5 * sum(mI_n - mP)
		psi_0 <- as.vector(B + 0.5 * (t(1 - vp) %*% (vy - mX %*% vm_beta)^2 + tr(mS_beta %*% (t(mX) %*% (mI_n - mP) %*% mX))))
	}

	return(list(mS_alpha=mS_alpha, mS_beta=mS_beta, vm_alpha=vm_alpha,
							vm_beta=vm_beta, veta=veta, vp=vp, mP=mP, r_1=r_1,
							psi_1=psi_1, r_0=r_0, psi_0=psi_0))
}

main <- function()
{
	data <- generate()
	result <- vb(data)
	print(result)
}
main()