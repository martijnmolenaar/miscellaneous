options(stringsAsFactors = FALSE)
library(ggplot2)
library(dplyr)
#library(spdep)
library(limSolve)
library(reshape2)
library(gpclib)


### functions
rotate_xy <- function(xy, angle, center = c(0, 0)) {
  xRot = center[1] + cos(angle) * (xy[, 1] - center[1]) - sin(angle) * (xy[, 2] - center[2])
  yRot = center[2] + sin(angle) * (xy[, 1] - center[1]) + cos(angle) * (xy[, 2] - center[2])
  cbind(xRot, yRot)
}


convex.poly <- function(nSides, area)
{
  radius <- sqrt((2 * area) / (nSides * sin((2 * pi) / nSides)))
  angle <- (2 * pi) / nSides
  
  radii <- rnorm(nSides, radius, radius / 10)
  angles <- rnorm(nSides, angle, angle / 10) * 1:nSides
  angles <- sort(angles)
  
  points <- list(x = NULL, y = NULL)
  points$x <- cos(angles) * radii
  points$y <- sin(angles) * radii
  
  m <- matrix(unlist(points), ncol = 2)
  m <- rbind(m, m[1, ])
  current.area <-
    0.5 * (sum(m[1:nSides, 1] * m[2:(nSides + 1), 2]) - sum(m[1:nSides, 2] *
                                                              m[2:(nSides + 1), 1]))
  
  points$x <- points$x * sqrt(area / current.area)
  points$y <- points$y * sqrt(area / current.area)
  
  return (points)
}

regular.poly <- function(nSides, area)
{
  radius <- sqrt((2 * area) / (nSides * sin((2 * pi) / nSides)))
  
  points <- list(x = NULL, y = NULL)
  angles <- (2 * pi) / nSides * 1:nSides
  
  points$x <- cos(angles) * radius
  points$y <- sin(angles) * radius
  
  return (points)
  
}


### generate 9 random cells
coordinates <- rep(seq(-2, 2, length.out = 3), 3)
x <- matrix(coordinates, ncol = 3, byrow = F)
y <- matrix(coordinates, ncol = 3, byrow = T)

cells <-
  sapply(1:9, function(i) {
    xy <- c(x[i], y[i])
    cell <- data.frame(cell = paste("cell#", i, sep = ""),
                       convex.poly(
                         nSides = round(runif(
                           n = 1, min = 15, max = 30
                         ), 0),
                         area = rnorm(n = 1, mean = 2.5, sd = .3)
                       ))
    cell[2:3] <- t(apply(cell[2:3], 1, function(row) {
      row + xy
    }))
    
    return(cell)
  }, simplify = FALSE) %>% do.call("rbind", .)

cells$cell <- factor(cells$cell, levels = paste("cell#", 1:9, sep = ""))


##### generate 16 ablations

coordinates <- rep(seq(-2, 2, length.out = 4), 4)
x <- as.vector(matrix(coordinates, ncol = 4, byrow = F))
y <- as.vector(matrix(coordinates, ncol = 4, byrow = T))

x <- rotate_xy(cbind(x, y), angle = 0.3)[, 1]
y <- rotate_xy(cbind(x, y), angle = 0.3)[, 2]

ablations <-
  sapply(1:16, function(i) {
    xy <- c(x[i], y[i])
    cell <- data.frame(ablation = paste("ablation#", i, sep = ""),
                       regular.poly(nSides = 40, area = .8))
    cell[2:3] <- t(apply(cell[2:3], 1, function(row) {
      row + xy
    }))
    
    return(cell)
  }, simplify = FALSE) %>% do.call("rbind", .)

ablations$ablation <-
  factor(ablations$ablation, levels = paste("ablation#", 1:16, sep = ""))

### plot simulation

ggplot() +
  coord_fixed() +
  geom_polygon(
    data = cells,
    aes(x = x, y = y, fill = cell),
    color = "black",
    size = 1
  ) +
  geom_polygon(
    data = ablations,
    aes(x = x, y = y, color = ablation),
    fill = "gray",
    alpha = .7,
    size = 2
  ) +
  theme_void()+
  theme(legend.direction = "vertical", 
        legend.position = "right",
        legend.box = "horizontal"
  )

#### assign random metabolite X signals

cell_signals <- data.frame(cell = paste("cell#", 1:9, sep = ""),
                           signal = runif(9, min = 0, max = 10))

cell_signals

cells$signal <-
  cell_signals$signal[match(cells$cell, cell_signals$cell)]

### plot signals

ggplot() +
  coord_fixed() +
  scale_fill_gradient(low = "gray", high = "red") +
  geom_polygon(
    data = cells,
    aes(
      x = x,
      y = y,
      group = cell,
      fill = signal
    ),
    color = "black",
    size = 1
  ) +
  geom_polygon(
    data = ablations,
    aes(x = x, y = y, group = ablation),
    color = "white",
    fill = "black",
    alpha = .9,
    size = 2,
    show.legend = F
  ) +
  theme_void()+
  theme(panel.background = element_rect(fill = "gray"))



###### calculate expected signal in the ablations

## first, formalize polys to measure them
cells_polys <-
  cells %>% group_split(cell) %>% sapply(function(cell) {
    shape <- cbind(cell$x, cell$y)
    shape <- shape[chull(shape),]
    shape <- as(shape, "gpc.poly")
    shape
  })

ablations_polys <-
  ablations %>% group_split(ablation) %>% sapply(function(ablation) {
    shape <- cbind(ablation$x, ablation$y)
    shape <- shape[chull(shape),]
    shape <- as(shape, "gpc.poly")
    shape
  })

## calculate area% of cell_i in abl_i
overlap_matrix <-
  sapply(1:length(ablations_polys), function(abl_i) {
    sapply(1:length(cells_polys), function(cell_i) {
      area.poly(intersect(ablations_polys[[abl_i]], cells_polys[[cell_i]]))
    }) / area.poly(ablations_polys[[abl_i]])
  })

overlap_matrix <- data.frame(overlap_matrix)
rownames(overlap_matrix) <- paste("cell#", 1:9, sep = "")
colnames(overlap_matrix) <- paste("ablation#", 1:16, sep = "")

## matrix looks like:
head(overlap_matrix)


## now calculate anticipated signals in ablations 
ablations_signals <-
  data.frame(
    ablation = paste("ablation#", 1:16, sep = ""),
    signal = apply(overlap_matrix, 2, function(abl_i) {
      sum(abl_i * cell_signals$signal)
    })
  )

ablations_signals

ablations$signal <-
  ablations_signals$signal[match(ablations$ablation, ablations_signals$ablation)]


## plot signals to check, seems to be ok
ggplot() +
  coord_fixed() +
  scale_fill_gradient(low = "gray", high = "red") +
  geom_polygon(
    data = cells,
    aes(
      x = x,
      y = y,
      group = cell,
      fill = signal
    ),
    color = "black",
    size = 1
  ) +
  geom_polygon(
    data = ablations,
    aes(
      x = x,
      y = y,
      group = ablation,
      fill = signal
    ),
    color = "white",
    size = 2,
    show.legend = F
  ) +
  theme_void()+
  theme(panel.background = element_rect(fill = "gray"))


### now calculate back with SpaceM method 
### the metabolite concentrations per cell
cell_signals$SpaceM <-
  sapply(1:9, function(cell_i) {
    top_bottom <- rowSums(sapply(1:16, function(abl_i) {
      top <-
        ablations_signals$signal[abl_i] * overlap_matrix[cell_i, abl_i] * 
        area.poly(intersect(ablations_polys[[abl_i]], cells_polys[[cell_i]]))
      bottom <-
        area.poly(intersect(ablations_polys[[abl_i]], cells_polys[[cell_i]]))
      c(top, bottom)
    }))
    top_bottom[1] / top_bottom[2]
  })

cell_signals

cell_signals_long <-
  melt(
    cell_signals,
    id.vars = 1:2 ,
    variable.name = "method",
    value.name = "back_estimation"
  )

## plot comparison
ggplot(data = cell_signals_long,
       aes(
         x = signal,
         y = back_estimation,
         color = cell,
         shape = method
       )) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  coord_fixed() +
  geom_point(size = 3) +
  geom_smooth(
    aes(linetype = method),
    method = "lm",
    se = F,
    color = "gray"
  ) +
  theme_classic()

### in fact, it is a linear inverse problem, with 16 
### knowns (ablations) and 10 unknowns (9 cells + background)
### A*x = B

B <- ablations_signals$signal
B

A <- overlap_matrix
A <- t(as.matrix(rbind(A, 1 - colSums(A)))) ## last row is background
colnames(A)[10] <- "background"
A

### Solve system:
cell_signals$LIMSolve <- Solve(A, B)[1:9] ##10th element is background
cell_signals

cell_signals_long <-
  melt(
    cell_signals,
    id.vars = 1:2 ,
    variable.name = "method",
    value.name = "back_estimation"
  )

### plot comparisons, LIMSolve seems to be more exact
ggplot(data = cell_signals_long,
       aes(
         x = signal,
         y = back_estimation,
         color = cell,
         shape = method
       )) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  coord_fixed() +
  geom_point(size = 3) +
  geom_smooth(
    aes(linetype = method),
    method = "lm",
    se = F,
    color = "gray"
  ) +
  theme_classic()
