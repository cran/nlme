### $Id: Tetracycline1.R,v 1.1 2000/03/17 22:21:21 saikat Exp $
### Tetracycline pharmacokinetic data
### Hand and Crowder (1996), Table B.8, p. 198
"Tetracycline1" <-
  structure(list
  (conc = c(1.08, 1.99, 1.46, 1.21, 1.48, 2.5, 2.62, 
     1.95, 1.19, 2.1, 1.21, 0.96, 0.62, 0.88, 0.68, 0.48, 1.22, 1.91, 
     1.36, 0.9, 0.65, 1.52, 1.32, 0.95, 0.6, 1.1, 1.03, 0.61, 0.32, 
     2.12, 1.48, 1.09, 0.55, 1, 0.82, 0.52, 1.48, 0.9, 0.75, 0.44),
   Time = c(1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3,
     6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 2, 3, 6),
   Subject = structure(ordered(c(5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3,
     3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1,
     1, 1, 1, 1, 1, 1), levels=1:5),
     class = c("ordered", "factor"),
     .Label = c("5", "3", "2", "4", "1")),
   Formulation = structure(factor(c(1, 
     1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 
     2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2), levels=1:2),
     class = "factor",
     .Label = c("tetrachel", "tetracyn"))),
row.names = 1:40,
class = c("nfnGroupedData", "nfGroupedData", "groupedData", "data.frame"),
formula = conc ~ Time | Subject,
labels = list(x = "Time since drug administration",
  y = "Serum concentration of Tetracycline"),
units = list(x = "(hr)", y = "(IU/ml)"),
inner = ~Formulation)
