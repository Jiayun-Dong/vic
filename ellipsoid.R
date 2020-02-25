ellipsoid = function(cen, len, rot, NUM){
  # sample from p-sphere
  # first sample unit-length vector with different directions
  points = matrix(runif(length(cen) * NUM, -1, 1), nrow = NUM, ncol = length(cen))
  r = sqrt(rowSums(points^2))
  points = diag(1/r) %*% points
  
  # sample distance from origin
  ## we can change the distribution so that more points are far from the origin
  # r = runif(NUM, 0 ,1)
  r = rbeta(NUM, 2, 1)
  # hist(r)
  # plot(density(r))
  
  # get points from p-sphere
  points = diag(r) %*% points
  
  # stretch, rotate and shift to get ellipsoid
  points = points %*% diag(len)
  points = points %*% t(rot)
  points = sweep(points, 2, cen, "+")
  
  return(points)
}