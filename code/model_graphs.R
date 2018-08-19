# model_graphs.R
require(igraph)

# Univariate model
univariate_edges = matrix(c("y", "r", "r", "rho", "y", "x", "x", "lambda"), ncol = 2, 
    byrow = TRUE)
g = graph.edgelist(univariate_edges)
layout.reingold.tilford(g, root = c(1))
plot(g)
tkplot(g)

# Multivariate model
multivariate_edges = matrix(c("y", "r", "r", "rho", "y", "nu", "nu", "beta", "nu", 
    "u", "u", "sigma2_u"), ncol = 2, byrow = TRUE)
g2 = graph.edgelist(multivariate_edges)
# layout.kamada.kawai(g2, root='y')
layout.reingold.tilford(g2, root = c(1))
plot(g2)
tkplot(g2)

g3 = graph.tree(12, children = 2)
plot(g3, layout = layout.reingold.tilford)

g3 = graph.tree(12, children = 2)
plot(g3, layout = layout.reingold.tilford)
