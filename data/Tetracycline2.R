### $Id: Tetracycline2.R,v 1.3 2000/07/03 18:22:45 bates Exp $
### Tetracycline pharmacokinetic data.
### Hand and Crowder (1996), Table B.8, p. 198
"Tetracycline2" <-
  structure(list
  (conc = c(1.2, 1.54, 1.28, 0.79, 1.28, 2.25, 1.95, 
     1.24, 0.96, 2.05, 1.65, 1.05, 0, 1.36, 1.24, 0.6, 1.89, 2.55, 
     2.35, 1.3, 1.42, 1.7, 1.42, 1.05, 0, 1.04, 0.94, 0.57, 0.64, 
     0.74, 0.5, 0, 0, 0.94, 0.96, 0.68, 0.92, 1.72, 1.65, 1.25),
   Time = c(1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3,
     6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6),  
   Subject = structure(ordered(c(4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3,
     3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
     2, 2, 2, 2, 2, 2), levels=1:5), class = c("ordered", "factor"),
     .Label = c("4", "5", "2", "1", "3")), 
   Formulation = structure(factor(c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1,
     1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1,
     1, 1, 1, 2, 2, 2, 2), levels=1:2), class = "factor",
     .Label = c("Berkmycin", "tetramycin"))), 
row.names = 1:40,
class = c("nfnGroupedData", "nfGroupedData", "groupedData", "data.frame"),
formula = conc ~ Time | Subject,
labels = list(x = "Time since drug administration",
  y = "Serum concentration of Tetracycline"),
units = list(x = "(hr)", y = "(IU/ml)"),
inner = ~Formulation)
