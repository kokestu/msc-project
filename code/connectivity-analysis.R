library(doParallel)

# Interested in how many patches we get with different coverage proportions.
n <- 20
trials <- 20
cov_range <- seq(0.05, 0.95, length.out=n)
coverage <- rep(cov_range, trials)
size <- 50

# Run in parallel for speed
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)
result <- foreach(i=1:(n*trials), .combine='rbind') %dopar%
{
  library(raster)
  library(sf)       # vector operations
  # Set up the field
  arena <- raster(xmn=0, xmx=size, ymn=0, ymx=size, res=1)
  crs(arena) <- NA
  values(arena) <- sample(
    c(0,1), size^2, replace=TRUE, prob=c(1-coverage[i], coverage[i])
  )
  polys <- rasterToPolygons(arena, dissolve=TRUE, fun=\(x) x==1)
  polys <- disaggregate(polys)
  polys <- as(polys, 'sf')
  polys$area <- st_area(polys)
  return(c(nrow(polys), mean(polys$area), coverage[i]))
}
stopCluster(cl)

# Looks like 0.7-0.8ish is where we start having approx a single patch
grouping <- factor(result[,3])
n_patches <- tapply(result[,1], grouping, mean)
avg_patch_size <- tapply(result[,2], grouping, mean)

pdf('../results/connectivity.pdf', height = 3.5)
# par(mfrow=c(2,1))
x <- result[,3]  # cov_range
y <- result[,1]  # n_patches
plot(
  x,y, xlim=c(0,1), ylim=c(0,max(y)),
  xlab='Coverage', ylab='Number of patches'
)
abline(v=0.9, col='red',lty=2)
# x <- cov_range
# y <- (avg_patch_size)
# plot(x,y, xlim=c(0,1), ylim=c(0,max(y)), xlab='Coverage', ylab='Log(average patch area)')
# lines(seq(min(x),max(x),length.out=100), seq(min(y),max(y),length.out=100),lty=2)
# abline(v=0.9, col='red',lty=2)
dev.off()
# x <- n_patches
# y <- sqrt(avg_patch_size)
# plot(
#   x,y, xlim=c(0,max(x)), ylim=c(0,max(y)),
#   xlab='Number of patches', ylab='Sqrt patch area'
# )
