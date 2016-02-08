if (! require(igraph)) {
  install.packages('igraph')
  library(igraph)
}
if (! require(assertthat)) {
  install.packages('assertthat')
  library(assertthat)
}
if (! require(Matrix)) {
  install.packages('Matrix')
  library(Matrix)
}
if (! require(rARPACK)) {
  install.packages('rARPACK')
  library(rARPACK)
}

load('network.RData')

# Arguments
#   n : sample size
#   p : probability of non-zero interaction
#   k : mean interactions, given non-zero
generate.messages = function(n, p, k)
{
  (runif(n) < p) * (1 + rgeom(n, 1 / (k+1)))
}

# Arguments
#   n   : numeric, size of each friend group
#   p_f : probability of zero interaction between friends
#   p_a : probability of zero interaction between acquaintances
#   k_f : mean number of non-zero interactions between friends
#   k_a : mean number of non-zero interactions betweena acquaintances
# Return
#   $A  : adjacency matrix
#   $S  : list, for each friend group, the people in it
generate.groups = function(n, p_f, p_a, k_f, k_a)
{
  pop = 1:sum(n)
  S = list()
  f = list()
  A = matrix(0, nrow = sum(n), ncol = sum(n))
  for (i in 1:length(n))
  {
    S[[i]] = sort(sample(pop, n[i]))
    pop = setdiff(pop, S[[i]])
    
    for (j in 1:i)
    {
      if (i == j)
      { # friends
        x = generate.messages(n[i] * (n[i] - 1) / 2, p_f, k_f)
        tmp = matrix(0, nrow = n[i], ncol = n[i])
        tmp[row(tmp) < col(tmp)] = x
        tmp = tmp + t(tmp)
      }
      else
      { # acquaintances
        x = generate.messages(n[i] * n[j], p_a, k_a)
        tmp = matrix(x, nrow = n[i], ncol = n[j])
      }
      
      A[S[[i]], S[[j]]] = tmp
      A[S[[j]], S[[i]]] = t(tmp)
      
      if (! all(A == t(A)))
        stop('A not symmetric: ', i, ' ', j)
    }
  }
  list(A = A, S = S)
}

# Arguments
#   A  : adjacency matrix
#   S  : list, for each friend group, the people in it
#   detailed : logical(1), if TRUE, returns a matrix of cut weights
#              between pairs of friend groups
cut.weight = function(A, S, detailed = F)
{
  assert_that(all(A == t(A)))

  S = .set.preprocess(S, nrow(A))
    
  cut = 0
  if (detailed)
    cut.mat = matrix(0, nrow = length(S), ncol = length(S))
  for (i in 2:length(S))
  {
    for (j in 1:(i-1))
    {
      cut = cut + sum(A[S[[i]], S[[j]]])
      if (detailed)
        cut.mat[i,j] = cut.mat[j,i] = sum(A[S[[i]], S[[j]]])
    }
  }
  
  if (detailed) return(cut.mat)
  cut
}

# Arguments
#   A       : adjacency matrix
#   S.color : (optional) list, for each friend group, the people in it
#   S.group : (optional) list, for each friend group, the people in it
plot.network = function(A, S.color = NULL, S.group = NULL)
{
  if (nrow(A) > 200) stop('Graph too large to plot')
  assert_that(all(A == t(A)))
  g = graph.adjacency(A, 'undirected', weighted = T)
  n = nrow(A)
  
  if (is.null(S.color) & is.null(S.group))
    plot.igraph(g, edge.width = log(E(g)$weight + 1))
  else if (is.null(S.color))
  {
    S.group = .set.preprocess(S.group, n)
    plot.igraph(g, edge.width = log(E(g)$weight + 1),mark.groups = S.group)
  }
  else
  {
    S.color = .set.preprocess(S.color, n)
    palette = rainbow(length(S.color))
    colors = rep('white', times = n)
    for (s in 1:length(S.color))
      colors[S.color[[s]]] = palette[s]
    
    if (is.null(S.group))
      plot.igraph(g, edge.width = log(E(g)$weight + 1), vertex.color = colors)
    else
      plot.igraph(g, edge.width = log(E(g)$weight + 1), mark.groups = S.group,
                  vertex.color = colors)
  }
}

# Arguments
#   v         : numeric, eigenvector of Laplacian
#   quantiles : (optional) numeric, in (0, 1)
# Return
#   list, for each friend group, the people in it
partition.split = function(v, quantiles = NULL)
{
  if (is.null(quantiles))
    return(.set.preprocess(list(which(v > 0)), length(v)))
  quantiles = c(0, round(length(v) * quantiles), length(v))
  v = order(v)
  S = list()
  for (i in 2:length(quantiles))
    S[[i-1]] = sort(v[(quantiles[i-1]+1):quantiles[i]])
  S
}

rand.index = function(S1, S2, n)
{
  S1 = .set.preprocess(S1, n)
  S2 = .set.preprocess(S2, n)

  M1 = .sameset.matrix(S1, n)
  M2 = .sameset.matrix(S2, n)

  M = M1 == M2
  sum(M[upper.tri(M)]) / choose(n, 2)
}

rand.index.network = function(S, network)
{
  if      (network == 1) M = .network$id1
  else if (network == 2) M = .network$id2
  else stop('network must = 1 or 2')
  n = nrow(M)
  S = .set.preprocess(S, n)
  M = M == .sameset.matrix(S, n)
  sum(M[upper.tri(M)]) / choose(n, 2)
}


.sameset.matrix = function(S, n)
{
  S = .set.preprocess(S, n)
  if (max(sapply(S, max)) != n)
    stop('S does match graph of size n = ', n)
  M = matrix(0, nrow = n, ncol = n)
  for (i in 1:length(S))
    M[S[[i]], S[[i]]] = 1
  M
}

.set.preprocess = function(S, n)
{
  assert_that(all(unlist(S, use.names = F) %in% 1:n))
  assert_that(! any(duplicated(unlist(S, use.names = F))))
  
#  if (sum(sapply(S, length)) < n)
#    S[[length(S) + 1]] = setdiff(1:n, unlist(S, use.names = F))
  
  S
}
