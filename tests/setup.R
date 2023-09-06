library(SAScii)
library(readr)

sas_url <-
	"https://www.cdc.gov/healthyyouth/data/yrbs/sadc_2019/2019-SADC-SAS-Input-Program.sas"

dat_url <-
	"https://www.cdc.gov/healthyyouth/data/yrbs/sadc_2019/sadc_2019_national.dat"
	
sas_positions <-
	parse.SAScii( sas_url )

sas_positions[ , 'varname' ] <-
	tolower( sas_positions[ , 'varname' ] )
	
variables_to_keep <-
	c( "sex" , "grade" , "race4" , "q30" , "year" , "psu" , "stratum" , "weight" )
	
sas_positions[ , 'column_types' ] <-
	ifelse( !( sas_positions[ , 'varname' ] %in% variables_to_keep ) , "_" ,
		ifelse( sas_positions[ , 'char' ] , "c" , "d" ) )

yrbss_tbl <-
	read_fwf(
		dat_url ,
		fwf_widths( 
			abs( sas_positions[ , 'width' ] ) , 
			col_names = sas_positions[ , 'varname' ] 
		) ,
		col_types = paste0( sas_positions[ , 'column_types' ] , collapse = "" ) ,
		na = c( "" , "." )
	)
	
yrbss_df <- data.frame( yrbss_tbl )
yrbss_df <- subset( yrbss_df , year %in% seq( 1991 , 2011 , 2 ) )

yrbss_df[ , 'ever_smoked' ] <-
	as.numeric( yrbss_df[ , 'q30' ] == 1 )

yrbss_df[ , 'q30' ] <- NULL
yrbss_df[ , 'sex' ] <- relevel( factor( yrbss_df[ , 'sex' ] ) , ref = "2" )

for ( i in c( 'race4' , 'grade' ) ){
	yrbss_df[ , i ] <- relevel( factor( yrbss_df[ , i ] ) , ref = "1" )
}
distinct_years_available <- length( seq( 1991 , 2011 , 2 ) )

# store the linear polynomials
c11l <- contr.poly( distinct_years_available )[ , ".L" ]

# store the quadratic polynomials
c11q <- contr.poly( distinct_years_available )[ , ".Q" ]

# store the cubic polynomials
c11c <- contr.poly( distinct_years_available )[ , ".C" ]
# year^1 term (linear)
yrbss_df[ , "t11l" ] <- c11l[ match( yrbss_df[ , "year" ] , seq( 1991 , 2011 , 2 ) ) ]

# year^2 term (quadratic)
yrbss_df[ , "t11q" ] <- c11q[ match( yrbss_df[ , "year" ] , seq( 1991 , 2011 , 2 ) ) ]

# year^3 term (cubic)
yrbss_df[ , "t11c" ] <- c11c[ match( yrbss_df[ , "year" ] , seq( 1991 , 2011 , 2 ) ) ]
options( survey.lonely.psu = "adjust" )

library(survey)

des <- 
	svydesign(
		id = ~psu , 
		strata = ~interaction( stratum , year ) ,
		data = yrbss_df , 
		weights = ~weight , 
		nest = TRUE
	)

prevalence_over_time <-
	svyby( 
		~ ever_smoked , 
		~ year , 
		des , 
		svymean , 
		na.rm = TRUE 
	)

# confirm prevalence rates match published estimates
# of high school students that ever smoked
stopifnot(
	all.equal( 
		round( coef( prevalence_over_time ) , 3 ) , 
		c( .701 , .695 , .713 , .702 , .704 , .639 , .584 , .543 , .503 , .463 , .447 ) ,
		check.attributes = FALSE
	)
)
linyear <- 
	svyglm(
		ever_smoked ~ sex + race4 + grade + t11l , 
		design = des , 
		family = quasibinomial
	)

summary( linyear )
quadyear <-
	svyglm(
		ever_smoked ~ sex + race4 + grade + t11l + t11q , 
		design = des , 
		family = quasibinomial 
	)

summary( quadyear )
cubyear <-
	svyglm(
		ever_smoked ~ sex + race4 + grade + t11l + t11q + t11c , 
		design = des , 
		family = quasibinomial 
	)
	
summary( cubyear )
marginals <- 
	svyglm(
		formula = ever_smoked ~ sex + race4 + grade ,
		design = des , 
		family = quasibinomial
	)
( means_for_joinpoint <- svypredmeans( marginals , ~factor( year ) ) )
# coerce the results to a data.frame object
means_for_joinpoint <- as.data.frame( means_for_joinpoint )

# extract the row names as the survey year
means_for_joinpoint[ , "year" ] <- as.numeric( rownames( means_for_joinpoint ) )

# must be sorted, just in case it's not already
means_for_joinpoint <- means_for_joinpoint[ order( means_for_joinpoint[ , "year" ] ) , ]
means_for_joinpoint[ , "wgt" ] <- with( means_for_joinpoint, ( mean / SE ) ^ 2 ) 
o <- lm( log( mean ) ~ year , weights = wgt , data = means_for_joinpoint )
library(segmented)

# find only one joinpoint
os <- segmented( o , ~year )

summary( os )
slope( os , APC = TRUE )
# find two joinpoints rather than only one
os2 <- segmented( o , ~year , npsi = 2 )

summary( os2 )

slope( os2 , APC = TRUE )
# calculate a five-timepoint linear contrast vector
c5l <- contr.poly( 5 )[ , 1 ]

# tack the five-timepoint linear contrast vectors onto the current survey design object
des <- update( des , t5l = c5l[ match( year , seq( 1991 , 1999 , 2 ) ) ] )

pre_91_99 <-
	svyglm(
		ever_smoked ~ sex + race4 + grade + t5l ,
		design = subset( des , year <= 1999 ) , 
		family = quasibinomial
	)

summary( pre_91_99 )

# confirm 1991-1999 trend coefficient matches published estimates
stopifnot( round( pre_91_99$coefficients['t5l'] , 5 ) == .03704 )
# calculate a seven-timepoint linear contrast vector
c7l <- contr.poly( 7 )[ , 1 ]

# tack the seven-timepoint linear contrast vectors onto the current survey design object
des <- update( des , t7l = c7l[ match( year , seq( 1999 , 2011 , 2 ) ) ] )

post_99_11 <-
	svyglm(
		ever_smoked ~ sex + race4 + grade + t7l ,
		design = subset( des , year >= 1999 ) , 
		family = quasibinomial
	)
	
summary( post_99_11 )

# confirm 1999-2011 trend coefficient matches published estimates
stopifnot( round( post_99_11$coefficients['t7l'] , 5 ) == -0.99165 )
