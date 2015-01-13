# make_accuracy_tables.R
require(stargazer)

# Table of accuracies for random intercepts model
load(file="accuracy_results_int.RData")
beta_table = cbind(var1$vbeta_accuracy, var2$vbeta_accuracy, var3$vbeta_accuracy, var4$vbeta_accuracy, var5$vbeta_accuracy)
u_table = cbind(var1$vu_accuracy, var2$vu_accuracy, var3$vu_accuracy, var4$vu_accuracy, var5$vu_accuracy)
sigma_table = cbind(var1$sigma2_u_accuracy, var2$sigma2_u_accuracy, var3$sigma2_u_accuracy, var4$sigma2_u_accuracy, var5$sigma2_u_accuracy)
rho_table = cbind(var1$rho_accuracy, var2$rho_accuracy, var3$rho_accuracy, var4$rho_accuracy, var5$rho_accuracy)
table1 = rbind(beta_table, u_table, sigma_table, rho_table)
rownames(table1) = c("vbeta_1", "vbeta_2", "vu_1", "vu_2", "vu_3", "vu_4", "vu_5",
                     "vu_6", "vu_7", "vu_8", "vu_9", "vu_10", "sigma2_u",
                     "rho")
colnames(table1) = c("laplace", "gva", "gva2", "gva2 fast", "gva_nr")
stargazer(round(table1, 3))

# Table of accuracies for random slopes model
load(file="accuracy_results_slope.RData")
beta_table = cbind(var1$vbeta_accuracy, var2$vbeta_accuracy, var3$vbeta_accuracy, var4$vbeta_accuracy, var5$vbeta_accuracy)
u_table = cbind(var1$vu_accuracy, var2$vu_accuracy, var3$vu_accuracy, var4$vu_accuracy, var5$vu_accuracy)
sigma_table = cbind(var1$sigma2_u_accuracy, var2$sigma2_u_accuracy, var3$sigma2_u_accuracy, var4$sigma2_u_accuracy, var5$sigma2_u_accuracy)
rho_table = cbind(var1$rho_accuracy, var2$rho_accuracy, var3$rho_accuracy, var4$rho_accuracy, var5$rho_accuracy)
table2 = rbind(beta_table, u_table, sigma_table, rho_table)
round(table2, 3)
