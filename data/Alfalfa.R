### $Id: Alfalfa.R,v 1.3 2000/07/03 18:22:44 bates Exp $
### Yields of three varieties of Alfalfa (T/Acre) in 1994 following third
###  date of cutting in 1943
### Snedecor and Cochran (1980), Table 16.15.1, p. 327
Alfalfa <-
  structure(list(
  Variety = structure(factor(c(2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3), levels=1:3),
    class = "factor", .Label = c("Cossack", "Ladak", "Ranger")),
Date = structure(ordered(c(1, 2, 3, 4, 1, 
  2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 
  3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 
  4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 
  1, 2, 3, 4), levels=1:4), class = c("factor"), 
  .Label = c("None", "S1", "S20", "O7")), 
Block = structure(factor(c(1, 1, 1, 1, 2, 
  2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 1, 1, 
  1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 
  6, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 
  6, 6, 6, 6), levels=1:6),
  class = "factor",
  .Label = c("1", "2", "3", "4", "5", "6")),
Yield = c(2.17, 1.58, 2.29, 2.23, 1.88, 
  1.26, 1.6, 2.01, 1.62, 1.22, 1.67, 1.82, 2.34, 1.59, 1.91, 2.1, 
  1.58, 1.25, 1.39, 1.66, 1.66, 0.94, 1.12, 1.1, 2.33, 1.38, 1.86, 
  2.27, 2.01, 1.3, 1.7, 1.81, 1.7, 1.85, 1.81, 2.01, 1.78, 1.09, 
  1.54, 1.4, 1.42, 1.13, 1.67, 1.31, 1.35, 1.06, 0.88, 1.06, 1.75, 
  1.52, 1.55, 1.56, 1.95, 1.47, 1.61, 1.72, 2.13, 1.8, 1.82, 1.99, 
  1.78, 1.37, 1.56, 1.55, 1.31, 1.01, 1.23, 1.51, 1.3, 1.31, 1.13, 
  1.33)),
row.names = 1:72,
class = c("nmGroupedData", "groupedData", "data.frame"),
formula = Yield ~ Date | Block/Variety,
formulaList = list(Block = ~Block, Variety = ~Variety),
order.groups = list(Block = TRUE, Variety = TRUE),
FUN = function (x) max(x, na.rm = TRUE),
labels = list(y = "Yield in 1944 following third date of cutting in 1943"),
units = list(y = "(T/Acre)"))

