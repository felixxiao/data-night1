source('lib.R')
set.seed(1)

## Simulate a social network
# specify one group of 12 friends, and one group of 9
n = c(12, 9)

res = generate.groups(n, 0.8, 0.1, 20, 3)
A = res$A  # adjacency matrix
S = res$S  # which people are in each friend group

# Laplacian matrix
L = diag(colSums(A)) - A
# 1st non-trivial eigenvector
v = eigs(L, 2, 'SM')$vectors[,1]

# create grouping of friends
X = partition.split(v)

# vertices are grouped according to S, colored according to X
plot.network(A, S.color = X, S.group = S)

# check the ratio of non-zero interactions within predicted groups
tmp = A[X[[1]], X[[1]]]
mean(tmp[upper.tri(tmp)] != 0)

# check ratio of non-zero interactions between different predicted groups
mean(A[X[[1]], X[[2]]] != 0)

cut.weight(A, S)
cut.weight(A, X)

# try with 3 simulated groups
n = c(10, 13, 14)
