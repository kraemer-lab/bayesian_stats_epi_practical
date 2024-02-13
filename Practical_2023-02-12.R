install.packages("ggplot2")
install.packages("raster")
install.packages("sf")
install.packages("tiff")
install.packages("viridis")
install.packages("devtools")
install.packages("rstan")
install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
devtools::install_github("rmcelreath/rethinking@slim")

library(raster)
library(sf)
library(tiff)
library(viridis)
library(rstan)
library(rethinking)
library(ggplot2)

# set your own working directory

setwd('/Users/mkraemer/Dropbox/Practical_2024-02-12/data')

chol_map <- read_sf("buildings/deaths_by_bldg.shp")
pump_map <- read_sf("pumps/pumps.shp")
snow_map <- brick('SnowMap.tif')
os_map <- brick('OSMap_Grayscale.tif')
os_map2 <- brick('OSMap.tif')

order_index <- order(chol_map$deaths)
chol_map_ordered <- chol_map[order_index, ]

chol_map_ordered$deaths <- factor(chol_map_ordered$deaths, levels=unique(chol_map_ordered$deaths))

color_palette <- viridis::plasma(length(unique(chol_map_ordered$deaths)))

palette(color_palette)

png(filename="js_scatter_colormap.png", units="in", width=8, height=8, res=200)
plotRGB(snow_map, interpolate=TRUE, maxpixels=500000000)
points(chol_map_ordered$COORD_X, chol_map_ordered$COORD_Y, col=chol_map_ordered$deaths, pch=19)
legend("topright", legend=unique(chol_map_ordered$deaths), fill=color_palette, title="Deaths")
dev.off()

chol_map$pump_id <- as.character(chol_map$pumpID)

chol_map$pump_idx <- as.numeric(factor(chol_map$pump_id, 
                                       levels=unique(chol_map$pump_id)))
deaths <- chol_map$deaths
distance <- chol_map$distBSpump
dz <- (distance - mean(distance)) / sd(distance) #standarised distance to improve sampling

data <- list(
  y = deaths,#cholera death counts per building
  d = dz, #distance to Broad Street pump (standarised, z-scores)
  p = chol_map$pump_idx, #id index of pump respective to building location
  h = chol_map$ID #id of each measured building 
)

data 

prior_model1 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a + b[p],
    a ~ normal( 0, 1 ),
    b[p] ~ normal(0 , 1 )
  ), data=data )


prior <- extract.prior(prior_model1, n=1000, pars=c(y))

y_ppc <- link(prior_model1, post=prior, data=data)

y_mean <- colMeans(y_ppc)

y_lim = max(y_mean)


plt <- ggplot() +
  xlim(0, 30) +
  ylim(0, 2) +
  labs(x = "Deaths (y)", y = "Density", title = "Model 1 Prior Predictive Check") +
  theme_minimal()

for (i in 500:600) {
  iteration_data <- data.frame(value = y_ppc[i, ])
  plt <- plt + geom_density(data = iteration_data, aes(x = value, color = "Prior Predictive Samples"), alpha = 0.05)
}

mean_data <- data.frame(value = colMeans(y_ppc))
plt <- plt + geom_density(data = mean_data, aes(x = value, color = "Prior Predictive Mean"), linetype = "dashed", size = 1)

obs_data <- data.frame(value = chol_map$deaths)
plt <- plt + geom_density(data = obs_data, aes(x = value, color = "Observed"), size = 1)

plt <- plt + scale_color_manual(name = "", 
                                breaks = c("Observed", "Prior Predictive Mean", "Prior Predictive Samples"), 
                                labels = c("Observed", "Prior Predictive Mean", "Prior Predictive Samples"),
                                values = c("black", "darkorange", "deepskyblue2"))  

plt <- plt + theme_bw() 

ggsave("model1_prior_predictives.png", plot = plt, width = 8, height = 8, dpi = 200, bg = NULL)


##sample model to get inference data (samples)
samp_model1 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a + b[p],
    a ~ normal( 0 , 1 ),
    b[p] ~ normal( 0 , 1 )
  ), data=data , chains=4 , log_lik=TRUE, iter=2000 )

#summary(samp_model1)

prec <- precis(samp_model1, depth=2, prob=0.9)
prec <- data.frame(prec)
write.csv(prec, "model1_summary.csv")

png(filename = "model1_trankplots.png", units = "in", width = 12, height = 5, res = 200)
par(oma = c(0, 0, 3, 0))
trankplot(samp_model1, n_cols=2, lwd=3)
mtext("Model 1 trace rank plots", side=3, line=1, at=0.5, cex = 1.5, outer=TRUE)
dev.off()
while (!is.null(dev.list()))  dev.off()

##extract posterior predictive and plot
post <- extract.samples(samp_model1, n=1000)
y_pred <- link(samp_model1, post=post, data=data)


plt <- ggplot() +
  xlim(0, 30) +
  ylim(0, 2) +
  labs(x = "Deaths (y)", y = "Density", title = "Model 1 Posterior Predictive Check") +
  theme_minimal()

for (i in 500:600) {
  iteration_data <- data.frame(value = y_pred[i, ])
  plt <- plt + geom_density(data = iteration_data, aes(x = value, color = "Posterior Predictive Samples"), alpha = 0.01)
}

mean_data <- data.frame(value = colMeans(y_pred))
plt <- plt + geom_density(data = mean_data, aes(x = value, color = "Posterior Predictive Mean"), linetype = "dashed", size = 1)

obs_data <- data.frame(value = chol_map$deaths)
plt <- plt + geom_density(data = obs_data, aes(x = value, color = "Observed"), size = 1)

plt <- plt + scale_color_manual(name = "", 
                                breaks = c("Observed", "Posterior Predictive Mean", "Posterior Predictive Samples"), 
                                labels = c("Observed", "Posterior Predictive Mean", "Posterior Predictive Samples"),
                                values = c("black", "purple", "limegreen"))  

plt <- plt + theme_bw() 

ggsave("model1_Posterior_predictives.png", plot = plt, width = 8, height = 8, dpi = 200, bg = NULL)


## plot posterior predictive mean on map
chol_map$preds <- round(colMeans(y_pred), 1)
order_index <- order(chol_map$preds)
chol_map_ordered <- chol_map[order_index, ]
chol_map_ordered$preds <- factor(chol_map_ordered$preds, levels=unique(chol_map_ordered$preds))
color_palette <- viridis::plasma(length(unique(chol_map_ordered$preds)))
palette(color_palette)

png(filename="model1_preds_map.png", units="in", width=8, height=8, res=200)
plotRGB(snow_map, interpolate=TRUE, maxpixels=500000000)
points(chol_map_ordered$COORD_X, chol_map_ordered$COORD_Y, col=chol_map_ordered$preds, pch=19)
legend("topright", legend=unique(chol_map_ordered$preds), fill=color_palette, title="Predicted Deaths")
dev.off()
while (!is.null(dev.list()))  dev.off()

##sample model to get inference data (samples)
samp_model1 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a + b[p],
    a ~ normal( 0 , 1 ),
    b[p] ~ normal( 0 , 1 )
  ), data=data , chains=4 , log_lik=TRUE, iter=2000 )

prec <- precis(samp_model1, depth=2, prob=0.9)
prec <- data.frame(prec)
write.csv(prec, "model1_summary.csv")
prec


####################### Model 2 ########################
#######################################################

##Compute and plot prior predictive checks
prior_model2 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a + b*d,
    a ~ normal( am , as ),
    am ~ normal(0, 1),
    as ~ exponential(1),
    b ~ normal( bm , bs ),
    bm ~ normal(0, 1), #hyperpriors
    bs ~ exponential(1) #hyperpriors
  ), data=data )

prior <- extract.prior(prior_model2, n=1000, pars=c(y))

y_ppc <- link(prior_model2, post=prior, data=data)

plt <- ggplot() +
  xlim(0, 30) +
  ylim(0, 2) +
  labs(x = "Deaths (y)", y = "Density", title = "Model 2 Prior Predictive Check") +
  theme_minimal()

for (i in 500:600) {
  iteration_data <- data.frame(value = y_ppc[i, ])
  plt <- plt + geom_density(data = iteration_data, aes(x = value, color = "Prior Predictive Samples"), alpha = 0.05)
}

mean_data <- data.frame(value = colMeans(y_ppc))
plt <- plt + geom_density(data = mean_data, aes(x = value, color = "Prior Predictive Mean"), linetype = "dashed", size = 1)

obs_data <- data.frame(value = chol_map$deaths)
plt <- plt + geom_density(data = obs_data, aes(x = value, color = "Observed"), size = 1)

plt <- plt + scale_color_manual(name = "", 
                                breaks = c("Observed", "Prior Predictive Mean", "Prior Predictive Samples"), 
                                labels = c("Observed", "Prior Predictive Mean", "Prior Predictive Samples"),
                                values = c("black", "darkorange", "deepskyblue2"))  

plt <- plt + theme_bw() 

ggsave("model2_prior_predictives.png", plot = plt, width = 8, height = 8, dpi = 200, bg = NULL)

##sample model to get inference data (samples)
samp_model2 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a + b*d,
    a ~ normal( am , as ),
    am ~ normal(0, 1),
    as ~ exponential(1),
    b ~ normal( bm , bs ),
    bm ~ normal(0, 1), #hyperpriors
    bs ~ exponential(1) #hyperpriors
  ), data=data , chains=4 , log_lik=TRUE, iter=2000 )

prec <- precis(samp_model2, depth=2, prob=0.9)
prec <- data.frame(prec)
write.csv(prec, "model2_summary.csv")

extract <- rstan::extract #to avoid conflict with other packages

png(filename = "model2_trankplots.png", units = "in", width = 12, height = 5, res = 200)
par(oma = c(0, 0, 3, 0))
trankplot(samp_model2, n_cols=2, lwd=3)
mtext("Model 2 trace rank plots", side=3, line=1, at=0.5, cex = 1.5, outer=TRUE)
dev.off()
while (!is.null(dev.list()))  dev.off()


##extract posterior predictive and plot
post <- extract.samples(samp_model2, n=1000)
y_pred <- link(samp_model2, post=post, data=data)

plt <- ggplot() +
  xlim(0, 30) +
  ylim(0, 2) +
  labs(x = "Deaths (y)", y = "Density", title = "Model 2 Posterior Predictive Check") +
  theme_minimal()

for (i in 500:600) {
  iteration_data <- data.frame(value = y_pred[i, ])
  plt <- plt + geom_density(data = iteration_data, aes(x = value, color = "Posterior Predictive Samples"), alpha = 0.01)
}

mean_data <- data.frame(value = colMeans(y_pred))
plt <- plt + geom_density(data = mean_data, aes(x = value, color = "Posterior Predictive Mean"), linetype = "dashed", size = 1)

obs_data <- data.frame(value = chol_map$deaths)
plt <- plt + geom_density(data = obs_data, aes(x = value, color = "Observed"), size = 1)

plt <- plt + scale_color_manual(name = "", 
                                breaks = c("Observed", "Posterior Predictive Mean", "Posterior Predictive Samples"), 
                                labels = c("Observed", "Posterior Predictive Mean", "Posterior Predictive Samples"),
                                values = c("black", "purple", "limegreen"))  

plt <- plt + theme_bw() 

ggsave("model2_Posterior_predictives.png", plot = plt, width = 8, height = 8, dpi = 200, bg = NULL)

## plot posterior predictive mean on map
post <- extract.samples(samp_model2, n=1000)
y_pred <- link(samp_model2, post=post, data=data)
chol_map$preds <- round(colMeans(y_pred), 1)
order_index <- order(chol_map$preds)
chol_map_ordered <- chol_map[order_index, ]
chol_map_ordered$preds <- factor(chol_map_ordered$preds, levels=unique(chol_map_ordered$preds))
color_palette <- viridis::plasma(length(unique(chol_map_ordered$preds)))
palette(color_palette)

png(filename="model2_preds_map.png", units="in", width=8, height=8, res=200)
plotRGB(snow_map, interpolate=TRUE, maxpixels=500000000)
points(chol_map_ordered$COORD_X, chol_map_ordered$COORD_Y, col=chol_map_ordered$preds, pch=19)
legend("topright", legend=unique(chol_map_ordered$preds), fill=color_palette, title="Predicted Deaths")
dev.off()
while (!is.null(dev.list()))  dev.off()

####################### Model 3 ########################
#######################################################

##Compute and plot prior predictive checks
prior_model3 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a[h] + b[h]*d,
    a[h] ~ normal( am , as ),
    am ~ normal(0, 1),
    as ~ exponential(1),
    b[h] ~ normal( bm , bs ),
    bm ~ normal(0, 1), #hyperpriors
    bs ~ exponential(1) #hyperpriors
  ), data=data )


prior <- extract.prior(prior_model3, n=1000, pars=c(y))

y_ppc <- link(prior_model3, post=prior, data=data)

plt <- ggplot() +
  xlim(0, 30) +
  ylim(0, 2) +
  labs(x = "Deaths (y)", y = "Density", title = "Model 3 Prior Predictive Check") +
  theme_minimal()

for (i in 500:600) {
  iteration_data <- data.frame(value = y_ppc[i, ])
  plt <- plt + geom_density(data = iteration_data, aes(x = value, color = "Prior Predictive Samples"), alpha = 0.05)
}

mean_data <- data.frame(value = colMeans(y_ppc))
plt <- plt + geom_density(data = mean_data, aes(x = value, color = "Prior Predictive Mean"), linetype = "dashed", size = 1)

obs_data <- data.frame(value = chol_map$deaths)
plt <- plt + geom_density(data = obs_data, aes(x = value, color = "Observed"), size = 1)

plt <- plt + scale_color_manual(name = "", 
                                breaks = c("Observed", "Prior Predictive Mean", "Prior Predictive Samples"), 
                                labels = c("Observed", "Prior Predictive Mean", "Prior Predictive Samples"),
                                values = c("black", "darkorange", "deepskyblue2"))  

plt <- plt + theme_bw() 

ggsave("model3_prior_predictives.png", plot = plt, width = 8, height = 8, dpi = 200, bg = NULL)


##sample model to get inference data (samples)
samp_model3 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a[h] + b[h]*d,
    a[h] ~ normal( am , as ),
    am ~ normal(0, 1),
    as ~ exponential(1),
    b[h] ~ normal( bm , bs ),
    bm ~ normal(0, 1), 
    bs ~ exponential(1)
  ), data=data , chains=4 , log_lik=TRUE, iter=2000 )

prec <- precis(samp_model3, depth=2, prob=0.9)
prec <- data.frame(prec)
write.csv(prec, "model3_summary.csv")


##sample model to get inference data (samples)
samp_model3 <- ulam(
  alist(
    ## this is the actual model
    y ~ poisson( lambda ),
    log(lambda) <- a + b*d,
    a <- al + az[h]*as,
    az[h] ~ normal( 0 , 1 ),
    al ~ normal(0, 1),
    as ~ exponential(1),
    b <- bl + bz[h]*bs,
    bz[h] ~ normal( 0 , 1 ),
    bl ~ normal(0, 1), 
    bs ~ exponential(1)
  ), data=data , chains=4 , log_lik=TRUE, iter=2000 )

prec <- precis(samp_model3, depth=2, prob=0.9)
prec <- data.frame(prec)
write.csv(prec, "model3_summary.csv")

png(filename = "model3_trankplots.png", units = "in", width = 12, height = 5, res = 200)
par(oma = c(0, 0, 3, 0))
trankplot(samp_model3, n_cols=2, lwd=3)
mtext("Model 3 trace rank plots", side=3, line=1, at=0.5, cex = 1.5, outer=TRUE)
while (!is.null(dev.list()))  dev.off()

##extract posterior predictive and plot
post <- extract.samples(samp_model3, n=1000)
y_pred <- link(samp_model3, post=post, data=data)$lambda


plt <- ggplot() +
  xlim(0, 30) +
  ylim(0, 2) +
  labs(x = "Deaths (y)", y = "Density", title = "Model 3 Posterior Predictive Check") +
  theme_minimal()

for (i in 500:600) {
  iteration_data <- data.frame(value = y_pred[i, ])
  plt <- plt + geom_density(data = iteration_data, aes(x = value, color = "Posterior Predictive Samples"), alpha = 0.01)
}

mean_data <- data.frame(value = colMeans(y_pred))
plt <- plt + geom_density(data = mean_data, aes(x = value, color = "Posterior Predictive Mean"), linetype = "dashed", size = 1)

obs_data <- data.frame(value = chol_map$deaths)
plt <- plt + geom_density(data = obs_data, aes(x = value, color = "Observed"), size = 1)

plt <- plt + scale_color_manual(name = "", 
                                breaks = c("Observed", "Posterior Predictive Mean", "Posterior Predictive Samples"), 
                                labels = c("Observed", "Posterior Predictive Mean", "Posterior Predictive Samples"),
                                values = c("black", "purple", "limegreen"))  

plt <- plt + theme_bw() 

ggsave("model3_Posterior_predictives.png", plot = plt, width = 8, height = 8, dpi = 200, bg = NULL)

## plot posterior predictive mean on map
chol_map$preds <- round(colMeans(y_pred), 1)
order_index <- order(chol_map$preds)
chol_map_ordered <- chol_map[order_index, ]
chol_map_ordered$preds <- factor(chol_map_ordered$preds, levels=unique(chol_map_ordered$preds))
color_palette <- viridis::plasma(length(unique(chol_map_ordered$preds)))
palette(color_palette)

png(filename="model3_preds_map.png", units="in", width=8, height=8, res=200)
plotRGB(snow_map, interpolate=TRUE, maxpixels=500000000)
points(chol_map_ordered$COORD_X, chol_map_ordered$COORD_Y, col=chol_map_ordered$preds, pch=19)
legend("topright", legend=unique(chol_map_ordered$preds), fill=color_palette, title="Predicted Deaths")
dev.off()


