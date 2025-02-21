# Fonction de seuillage doux (soft-thresholding)
soft_thresholding <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}