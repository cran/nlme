### $Id: zzMethods.R,v 1.5.2.1 2003/03/25 15:23:54 bates Exp $
###
###   Miscellaneous methods that must be defined last in the library
###
### Copyright 1997-2001  Jose C. Pinheiro <jcp@research.bell-labs.com>,
###                      Douglas M. Bates <bates@stat.wisc.edu>
###
### This file is part of the nlme library for S and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

## Note that  require( nls )  has already happened ...

BIC.lme <- BIC.lmList <- BIC.gls <- BIC.lm

comparePred.lme <- comparePred.lmList <- comparePred.gls

getData.nlme <- getData.gnls

getData.lme <- getData.gls <- getData.nls

qqnorm.gls <- qqnorm.lm <- qqnorm.nls

plot.lme <- plot.nls

fitted.gnls <- fitted.gls

residuals.gnls <- residuals.gls

