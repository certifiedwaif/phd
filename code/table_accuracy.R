# make_accuracy_tables.R

load(file="accuracy_results_slope.RData")
beta_table = cbind(var1$vbeta_accuracy, var2$vbeta_accuracy, var3$vbeta_accuracy, var4$vbeta_accuracy, var5$vbeta_accuracy)
u_table = cbind(var1$vu_accuracy, var2$vu_accuracy, var3$vu_accuracy, var4$vu_accuracy, var5$vu_accuracy)
sigma_table = cbind(var1$sigma2_u_accuracy, var2$sigma2_u_accuracy, var3$sigma2_u_accuracy, var4$sigma2_u_accuracy, var5$sigma2_u_accuracy)
rho_table = cbind(var1$rho_accuracy, var2$rho_accuracy, var3$rho_accuracy, var4$rho_accuracy, var5$rho_accuracy)
table2 = rbind(beta_table, u_table, sigma_table, rho_table)
round(table2, 3)
