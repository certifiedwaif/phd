# regression_mixture.R

for (i in 1:ITERATIONS) {
	m_S_alpha = solve(r_1/psi_1)*t(mX)%*%mP%*%mX + sigmainv2_alpha * I)
	vm_alpha = m_S_alpha %*% t(mX) %*% mP %*% vy * r_1/psi_1
	m_S_beta = solve(r_0/psi_0)*t(mX)%*%mP%*%mX + sigmainv2_beta * I)
	vm_beta = m_S_beta %*% t(mX) %*% mP %*% vy * r_0/psi_0
	veta = logit(rho) + (-0.5*log(2*pi) - 0.5*log(psi_1) + 0.5*digamma(r_1) - 0.5*r_1/psi_1*((vy - t(mX)%*%vm_alpha)^2 + t(mX)%*%mS_alpha%*%mX))
	veta = veta - (-0.5*log(2*pi) - 0.5*log(psi_0) + 0.5*digamma(r_0) - 0.5*r_0/psi_0*((vy - t(mX)%*%vm_beta)^2 + t(mX)%*%mS_beta%*%mX))
	vp = expit(veta)
	mP = diag(vp)
	r_1 = A_1 + 0.5 * sum(vp)
	psi_1 = B + 0.5 * (t(vp) %*% (vy - mX%*%vm_alpha)^2 + tr(m_S_alpha %*% t(mX) %*% mP %*% mX))
	r_1 = A_2 + 0.5 * sum(vp)
	psi_1 = B + 0.5 * (t(vp) %*% (vy - mX%*%vm_beta)^2 + tr(m_S_beta %*% t(mX) %*% mP %*% mX))
}