### $Id: Meat.R,v 1.3 2000/07/03 18:22:45 bates Exp $
### Meat tenderness as a function of Storage period and Block
### Cochran and Cox (1957), section 11.51, p. 444
"Meat" <-
  structure(list
  (Storage = structure(factor(c(1, 2, 3, 4, 5, 6, 
     1, 3, 2, 5, 4, 6, 1, 4, 2, 6, 3, 5, 1, 5, 2, 4, 3, 6, 1, 6, 2, 
     3, 4, 5), levels=1:6), class = c("ordered", "factor"),
     .Label = c(" 0", " 1", " 2", " 4", " 9", "18")),
   score = c(7, 17, 26, 25, 33, 
     29, 17, 27, 23, 27, 29, 30, 10, 25, 26, 37, 24, 26, 25, 40, 25, 
     34, 34, 32, 11, 27, 24, 21, 26, 32),
   Block = structure(factor(c(3, 3, 3, 3, 3, 3, 1, 1, 1,
     1, 1, 1, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 2, 2, 2,
     2, 2, 2), levels=1:5),
     class = c("ordered", "factor"),
     .Label = c("II", "V", "I", "III", "IV")), 
   Pair = structure(ordered(c(7, 7, 8, 8, 9, 9, 1, 1, 2, 2, 3, 3, 10,
     10, 12, 12, 11, 11, 15, 15, 13, 13, 14, 14, 5, 5, 4, 4, 6, 6),
     levels=1:15), class = c("ordered", "factor"),
     .Label = c("II-1", "II-2", "II-3", "V-2", "V-1", "V-3", "I-1",
       "I-2", "I-3", "III-1", "III-3", "III-2", "IV-2", "IV-3",
       "IV-1"))),
row.names = 1:30,
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
formula = score ~ 1 | Pair,
outer = ~ Block,
labels = list(x = "Tenderness score", y = "Pair"),
inner = ~ Storage)
