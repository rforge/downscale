### Hypothetical species occupancy data
occ <- c(0.16, 0.36, 0.59, 0.86)
### cell areas (km2) for observed grain sizes of hypothetical species
areas <- c(40, 160, 640, 2560)

### Nachman model
(nachman <- downscale(occupancy = occ, area = areas, model = "Nachman"))

### Finite negative binomial model
(fnb <- downscale(occupancy = occ, area = areas, model = "FNB", 
                  extent = 256000))

### Thomas model
(thomas <- downscale(occupancy = occ, area = areas, model = "Thomas", 
                     extent = 256000, tolerance = 1e-5))

### grain sizes for prediction
fine.areas <- c(1, 2.5, 5, 10, 20, 40, 160, 640, 2560)

### model predictions
predict(nachman, newdata = fine.areas, plot = TRUE)
predict(fnb, newdata = fine.areas, plot = TRUE,
        extent = 256000)
predict(thomas, newdata = fine.areas, plot = TRUE, 
        extent = 256000, tolerance = 1e-5)

### ensemble modelling
ensemble.downscale(occupancy = occ, 
                   area = areas,
                   newdata = fine.areas,
                   extent = 256000,
                   tolerance_mod = 1e-5,
                   tolerance_pred = 1e-5,
                   model = "all")

