
///////////////////////////////////////////////////////////////////////////////////////////////////
// START WITH PULLING DATA FROM ONLINE SOURCES AND COMPILING IT INTO A SINGLE ANALYSABLE DATASET //
///////////////////////////////////////////////////////////////////////////////////////////////////

clear


// IMPORTANT: Start in the Folder that contains the Folders "Compiled Files" and "Results Files"
cd "Compiled Files"

// get global population size, density, mortality, and sex composition data
import delimited "https://raw.githubusercontent.com/cjm972/public-transit-CoV/main/CSV%20Data/WPP2022_Demographic_Indicators_Medium.csv"

keep if loctypename == "Country/Area"
rename tpopulation1july pop_size
rename popdensity pop_density
rename popsexratio pop_sexratio
rename medianagepop median_age

// we define a variable that indicates whether a country is part of the EEA (or in Switzerland's case part of the European single market)
gen EEA = 0

foreach i in "Austria" "Belgium" "Bulgaria" "Croatia" "Cyprus" "Czechia" "Denmark" "Estonia" "Finland" "France" "Germany" "Greece" "Hungary" "Ireland" "Italy" "Latvia" "Lithuania" "Luxembourg" "Malta" "Netherlands" "Poland" "Portugal" "Romania" "Slovakia" "Slovenia" "Spain" "Sweden" "Norway" "Liechtenstein" "Iceland" "Switzerland" {
	replace EEA = 1 if location=="`i'"
}
rename location country

// drop all observations not in the EEA + Switzerland group (EU+ group)
drop if EEA==0

// shift deaths by 1 year, 2 years, and 3 years to get deaths in 2019, 2020, 2021
gen deaths2019 = deaths[_n+3]
gen deaths2020 = deaths[_n+2]
gen deaths2021 = deaths[_n+1]

// drop all observations not from 2022
keep if time == 2022
rename time year


// drop all irrelevant variables
keep country year pop_size pop_density pop_sexratio median_age deaths2019 deaths2020 deaths2021

// save data as population dataset
save population.dta, replace
clear


// get global population age bucket data
import delimited "https://raw.githubusercontent.com/cjm972/public-transit-CoV/main/CSV%20Data/WPP2022_PopulationByAge5GroupSex_Percentage_Medium.csv"

keep if loctypename == "Country/Area"
rename poptotal pop_share

// we define a variable that indicates whether a country is part of the EEA (or in Switzerland's case part of the European single market)
gen EEA = 0

foreach i in "Austria" "Belgium" "Bulgaria" "Croatia" "Cyprus" "Czechia" "Denmark" "Estonia" "Finland" "France" "Germany" "Greece" "Hungary" "Ireland" "Italy" "Latvia" "Lithuania" "Luxembourg" "Malta" "Netherlands" "Poland" "Portugal" "Romania" "Slovakia" "Slovenia" "Spain" "Sweden" "Norway" "Liechtenstein" "Iceland" "Switzerland" {
	replace EEA = 1 if location=="`i'"
}
rename location country

// drop all observations not in the EEA + Switzerland group
drop if EEA==0

// drop all observations not from 2022
keep if time == 2022

// generate 6 20-year span age buckets
gen age_group = "none_"
replace age_group = "0to19" if agegrp == "0-4" | agegrp == "5-9" | agegrp == "10-14" | agegrp == "15-19"
replace age_group = "20to39" if agegrp == "20-24" | agegrp == "25-29" | agegrp == "30-34" | agegrp == "35-39"
replace age_group = "40to59" if agegrp == "40-44" | agegrp == "45-49" | agegrp == "50-54" | agegrp == "55-59"
replace age_group = "60to79" if agegrp == "60-64" | agegrp == "65-69" | agegrp == "70-74" | agegrp == "75-79"
replace age_group = "80to99" if agegrp == "80-84" | agegrp == "85-89" | agegrp == "90-94" | agegrp == "95-99"
replace age_group = "100plus" if agegrp == "100+"

// sum population share across the 6 age buckets and collapse the dataset to these 6 age buckets
collapse (sum) pop_share, by(country age_group)

// drop all irrelevant variables
keep country age_group pop_share

// just renaming pop_share into age_ for aesthetic purposes
rename pop_share age_

// reshape into correct format for merging later
reshape wide age_, i(country) j(age_group) string

// we merge the previous population data with this dataset to add the age data
merge 1:1 country using population.dta, nogenerate
save population.dta, replace
clear



// now load in GDP per capita 2021 data
import delimited using "https://raw.githubusercontent.com/cjm972/public-transit-CoV/main/CSV%20Data/Results.csv"
rename countryarea country
rename gdppercapitagdpusdollars gdppercap

// we use the most recent data that is fully available from the UN which is 2021 (the year before the time period our analysis investigates)
drop if year != 2021

// we define a variable that indicates whether a country is part of the EEA (or in Switzerland's case part of the European single market)
gen EEA = 0

foreach i in "Austria" "Belgium" "Bulgaria" "Croatia" "Cyprus" "Czechia" "Denmark" "Estonia" "Finland" "France" "Germany" "Greece" "Hungary" "Ireland" "Italy" "Latvia" "Lithuania" "Luxembourg" "Malta" "Netherlands" "Poland" "Portugal" "Romania" "Slovakia" "Slovenia" "Spain" "Sweden" "Norway" "Liechtenstein" "Iceland" "Switzerland" {
	replace EEA = 1 if country=="`i'"
}

// drop all observations not in the EEA + Switzerland group
drop if EEA==0
drop EEA

// drop irrelevant variables
keep country gdppercap

// merge into population.dta dataset (since the formats match and we keep most controls there)
merge 1:1 country using population.dta, nogenerate

// log transform gdp per capita
gen lgdp_cap = log(gdppercap)
drop gdppercap

save population.dta, replace
clear




// get global covid data
import delimited "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"

rename lat latitude
rename v4 longitude

// note that only Australia, China, and Canada are split up by province/state, some European countries have additional data for overseas territories but those don't concern us since they are not actually on the European continent.
// Thus, we can discard all data that contains province/state values
drop if provincestate != ""
drop provincestate

// we define a variable that indicates whether a country is part of the EEA (or in Switzerland's case part of the European single market)
gen EEA = 0

foreach i in "Austria" "Belgium" "Bulgaria" "Croatia" "Cyprus" "Czechia" "Denmark" "Estonia" "Finland" "France" "Germany" "Greece" "Hungary" "Ireland" "Italy" "Latvia" "Lithuania" "Luxembourg" "Malta" "Netherlands" "Poland" "Portugal" "Romania" "Slovakia" "Slovenia" "Spain" "Sweden" "Norway" "Liechtenstein" "Iceland" "Switzerland" {
	replace EEA = 1 if countryregion=="`i'"
}

// drop all observations not in the EEA + Switzerland group
drop if EEA==0
drop EEA latitude longitude


// transform data so as to have columns for each country and rows for each date
forval i = 1/31 {
	local label`i' = countryregion[`i']
}

forval j = 5/1047 {
	local dates`j' : var label v`j'
	disp("`dates`j''")
}

drop countryregion

xpose, clear

forval i = 1/31 {
	//replace v`i' = v`i' - v`i'[_n-1]
	rename v`i' `label`i''
	//local newname`i' = "covidcases_" + "`label`i''"
	local newname`i' = "covidcases" + "`i'"
	rename `label`i'' `newname`i''
}

gen report_date = ""
forval j = 5/1047 {
	replace report_date = "`dates`j''"  if _n==`j'-4
}

// generate a variable for date, year, and week
gen date = date(report_date, "MD20Y")
format date %tdDDmonCCYY
drop report_date
gen year = year(date)
gen week = week(date)
duplicates drop date, force

// reshape into correct format for merging
reshape long covidcases, i(date) j(country)

// make covid cases non-cumulative
by country (date), sort: gen diff_covidcases = covidcases - covidcases[_n-1]
drop covidcases
rename diff_covidcases covidcases
drop if year < 2022

// label countries correctly
tostring country, gen(country_string)
drop country
rename country_string country
forval i = 1/31 {
	disp("`label`i''")
	replace country = "`label`i''"  if country=="`i'"
}



// collapse to look at the sum of cases per week
collapse (mean) covidcases, by(country year week)

// lag outcome variables and drop the final weeks (where we don't have data)
gen cases_lagged = covidcases[_n+2]
drop if week > 46

// we merge the population data with this dataset
merge m:1 country year using population.dta, nogenerate

// normalize lagged cases by the population size
replace cases_lagged = cases_lagged/(pop_size*1000)

save global_covid.dta, replace
clear



// get government COVID-19 response stringency data
import delimited using "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/timeseries/government_response_index_avg.csv"

replace country_name = "Czechia" if country_name == "Czech Republic"
replace country_name = "Slovakia" if country_name == "Slovak Republic"

// note that only Australia, China, and Canada are split up by province/state, some European countries have additional data for overseas territories but those don't concern us since they are not actually on the European continent.
// Thus, we can discard all data that contains province/state values
drop if region_name != "NA"
drop region_name

// we define a variable that indicates whether a country is part of the EEA (or in Switzerland's case part of the European single market)
gen EEA = 0

foreach i in "Austria" "Belgium" "Bulgaria" "Croatia" "Cyprus" "Czechia" "Denmark" "Estonia" "Finland" "France" "Germany" "Greece" "Hungary" "Ireland" "Italy" "Latvia" "Lithuania" "Luxembourg" "Malta" "Netherlands" "Poland" "Portugal" "Romania" "Slovakia" "Slovenia" "Spain" "Sweden" "Norway" "Liechtenstein" "Iceland" "Switzerland" {
	replace EEA = 1 if country_name=="`i'"
}

// drop all observations not in the EEA + Switzerland group
drop if EEA==0
drop EEA

// drop all irrelevant variables
drop v1 jan2020 region_code jurisdiction country_code
forval i = 8/27 {
	drop v`i'
}

rename country_name countryregion

// rename some variables so that dates have consistent names
rename feb2020 v38
rename mar2020 v67
rename apr2020 v98
rename may2020 v128
rename jun2020 v159
rename jul2020 v189
rename aug2020 v220
rename sep2020 v251
rename oct2020 v281
rename nov2020 v312
rename dec2020 v342

rename jan2021 v373
rename feb2021 v404
rename mar2021 v432
rename apr2021 v463
rename may2021 v493
rename jun2021 v524
rename jul2021 v554
rename aug2021 v585
rename sep2021 v616
rename oct2021 v646
rename nov2021 v677
rename dec2021 v707

rename jan2022 v738
rename feb2022 v769
rename mar2022 v797
rename apr2022 v828
rename may2022 v858
rename jun2022 v889
rename jul2022 v919
rename aug2022 v950
rename sep2022 v981
rename oct2022 v1011
rename nov2022 v1042


// transform data so as to have columns for each country and rows for each date
forval i = 1/31 {
	local label`i' = countryregion[`i']
}

forval j = 28/1070 {
	local dates`j' : var label v`j'
	disp("`dates`j''")
}

drop countryregion

xpose, clear

forval i = 1/31 {
	//replace v`i' = v`i' - v`i'[_n-1]
	rename v`i' `label`i''
	//local newname`i' = "covidcases_" + "`label`i''"
	local newname`i' = "govstringency" + "`i'"
	rename `label`i'' `newname`i''
}

gen report_date = ""
forval j = 28/1070 {
	replace report_date = "`dates`j''"  if _n==`j'-27
}

// generate a variable for date, year, and week
gen date = date(report_date, "DMY")
format date %tdDDmonCCYY
drop report_date
gen year = year(date)
gen week = week(date)

// ensure correct time range is being looked at by dropping every year that is not 2022 and every week in 2022 beyond the 46th (time of Malta's policy)
drop if year != 2022
drop if week > 46

// reshape into correct format for merging
reshape long govstringency, i(date) j(country)

// label countries correctly
tostring country, gen(country_string)
drop country
rename country_string country
forval i = 1/31 {
	disp("`label`i''")
	replace country = "`label`i''"  if country=="`i'"
}

// collapse to look at the average government stringency per week
collapse (mean) govstringency, by(country year week)



// we merge the covid data with this dataset
merge 1:1 country year week using global_covid.dta, nogenerate

save global_covid.dta, replace



clear
// we now load in vaccination data
import delimited using "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv"
rename location country
rename total_vaccinations_per_hundred tot_vacc_rate
rename people_vaccinated_per_hundred people_vacc_rate
replace tot_vacc_rate = tot_vacc_rate/100
replace people_vacc_rate = people_vacc_rate/100

// only keep relevant variables
keep date country tot_vacc_rate  people_vacc_rate

// generate date variables
rename date rawdate
gen date = date(rawdate, "YMD")
format date %tdDDmonCCYY
drop rawdate
gen year = year(date)
gen week = week(date)

drop if year != 2022
drop if week > 46


// we define a variable that indicates whether a country is part of the EEA (or in Switzerland's case part of the European single market)
gen EEA = 0

foreach i in "Austria" "Belgium" "Bulgaria" "Croatia" "Cyprus" "Czechia" "Denmark" "Estonia" "Finland" "France" "Germany" "Greece" "Hungary" "Ireland" "Italy" "Latvia" "Lithuania" "Luxembourg" "Malta" "Netherlands" "Poland" "Portugal" "Romania" "Slovakia" "Slovenia" "Spain" "Sweden" "Norway" "Liechtenstein" "Iceland" "Switzerland" {
	replace EEA = 1 if country=="`i'"
}

// drop all observations not in the EEA + Switzerland group
drop if EEA==0
drop EEA

// collapse to look at the max vaccination rates per week
collapse (max) tot_vacc_rate people_vacc_rate, by(country year week)
drop if week > 46

// we merge the covid data with this dataset
merge 1:1 country year week using global_covid.dta, nogenerate

// NOTE: we will use total vaccination rates instead of total people vaccinated since we have less missing values for total vaccination rates

save global_covid.dta, replace
clear matrix
clear



//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//////////////////////////////////////////////
// START WITH ANALYSIS OF COMPILED DATASETS //
//////////////////////////////////////////////


// Using our compiled datasets, we now create a synthetic Germany composed of several untreated countries from the EU+ donor group and compare this synthetic untreated Germany to the real Germany to evaluate whether the treatment (9-Euro-Ticket over the summer) had any effect on the target (lagged Covid-19 cases). We construct our synthetic Germany using the synth package which allows us to specify the pre-treatment period on which to optimally match a combination of donor states to real Germany's actual lagged Covid-19 cases and the other control variables to optimally match across the entire time period investigated.

use global_covid.dta
cd ..
// switch file saving location
cd "Result Files"


// drop Luxembourg since Luxembourg had similar public transport subsidy to 9-Euro ticket and thus shouldn't be used as an untreated donor unit
drop if country == "Luxembourg"
// since Malta implements a transportation subsidy starting October 1st, 2022, we drop all observations past October 1st, 2022
drop if week > 39


// generate a numeric identifier for each country (necessary for the panel data setup we need to use in the synth package)
encode country, generate(country_id)

// Find the number by which Germany is coded in the dataset
gen temp_id = country_id if country == "Germany"
qui sum temp_id
local germany_id = `r(max)' 

// germany_id contains the number coding Germany
drop temp_id


// QUICK DIFFERENCES IN DIFFERENCES ANALYSIS TO COMPLEMENT SYNTHETIC CONTROL ANALYSIS

// Before we start the full synthetic control analysis, we do a quick Differences in Differences estimation and save it for comparison with our synthetic control analysis

// install outreg2 to save DiD results
ssc install outreg2

// generate a binary treated variable for DiD results for time period where 9 Euro ticket was implemented in Germany
gen treated = country == "Germany" & week >= 22 & week < 34

didregress (cases_lagged tot_vacc_rate govstringency age_0to19 age_100plus age_20to39 age_40to59 age_60to79 age_80to99 pop_density pop_sexratio lgdp_cap) (treated), group(country_id) time(week)

// save results in did_results_all.txt file in Results folder
outreg2 using did_results_all, text stats(coef se tstat pval) replace

// drop treated variable since we don't need it
drop treated


// plot lagged Covid-19 cases of Germany overlaid with the EU+ mean cases

// Note that week 22 was when treatment started (i.e., when the 9-Euro-Ticket was implemented) and week 34 was when it ended
// find average cases for entire donor pool and store them in avg_cases variable
bysort country week : egen avg_cases=mean(cases_lagged)
// seperate average cases into 31 seperate pieces equally long (but all containing the same EU+ mean lagged Covid-19 case numbers), this allows us to plot them on same plot as the cases by country
seperate avg_cases, by(country_id)
// seperate cases_lagged by country (useful for plotting lagged Covid-19 cases by country)
seperate cases_lagged, by(country_id)


// we plot the average amount of lagged cases for EU+ donor pool (avg_cases1, could pick avg_cases2 or any other avg_cases variable) and the actual lagged cases for Germany to compare
twoway line avg_cases1 week, lpattern(dash) lcolor(gray) || line cases_lagged`germany_id' week, ytitle("Lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) legend(order(1 "EU+ Mean" 2 "Germany") position(6) rows(1) region(lcolor(black))) lpattern(solid) lcolor(black) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(solid) lcolor(red))

// save plot of EU+ mean lagged Covid-19 cases overlaid on German Covid-19 cases
graph save eu_mean_all, replace
graph export eu_mean_all.eps, replace


// install synthetic control package and mat2txt package (will be useful for generating .txt files and .csv files for saving our results)
ssc install synth, replace all
ssc install mat2txt, replace all
gen order = 0


// set countries to units and weeks to time units (thus setting up the necessary panel data format for the synth package)
tsset country_id week


// Create synthetic Germany from EU+ donor pool
synth cases_lagged cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus, trunit(`germany_id') trperiod(22) margin(0.05) fig

// Uncomment to save graph of synthetic Germany generated by the synth function (we generate the same graph in several lines again just using the predicted synthetic cases directly)
graph save synthetic_germany_all_synth, replace
graph export synthetic_germany_all_synth.eps, replace


// Plot real Germany's lagged Covid-19 cases against those of the donor countries comprising synthetic Germany
twoway line cases_lagged`germany_id' week, lpattern(solid) lcolor(black) lwidth(thick) || line cases_lagged16 week, lpattern(dash) lcolor(blue) || line cases_lagged1 week, lpattern(dash) lcolor(red) || line cases_lagged9 week, lpattern(dash) lcolor(emerald) || line cases_lagged10 week, lpattern(dash) lcolor(gold) || line cases_lagged20 week, lpattern(dash) lcolor(brown) ytitle("Lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) legend(order(1 "Germany" 2 "Italy" 3 "Austria" 4 "Finland" 5 "France" 6 "Malta")  position(6) rows(3) region(lcolor(black))) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(solid) lcolor(red))

// save graph of synthetic components lagged Covid-19 cases versus Germany's lagged Covid-19 cases
graph save synthetic_components_all, replace
graph export synthetic_components_all.eps, replace


// Create a balance table (showing the mean values across our time period for each of the control and target variables for both synthetic and real Germany and the EU+ donor pool)
// loop through each variable used to create synthetic Germany (both controls and targets) and find their mean for entire EU+ group
foreach i in cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus {
	sum `i'
	local mean`i' = r(mean)
}

// build a matrix containing the average values for each variable for the donor pool
matrix eu_avg_mat = (`meancases_lagged'\ `meantot_vacc_rate'\ `meanpop_sexratio'\ `meangovstringency'\ `meanpop_density'\ `meanlgdp_cap'\ `meanage_0to19'\ `meanage_20to39'\ `meanage_40to59'\ `meanage_60to79'\ `meanage_80to99'\ `meanage_100plus')
// name the matrix column
matrix colnames eu_avg_mat = "EU+ Mean"
// save the balance table generated by the synth function as a matrix containing the variable values for synthetic and real Germany
matrix balance1 = e(X_balance)
// combine the matrices containing variable means for synthetic and real Germany and variable means for the EU+ group
matrix balance = balance1, eu_avg_mat
// save the final balance matrix with all 3 columns as a .txt file (can be thus converted into .csv very easily)
mat2txt, matrix(balance) saving(balance_all.txt) replace


// Save the weights table generated by the synth function as a matrix containing the weights of each individual donor to synthetic Germany
matrix country_weights = e(W_weights)
// save the weights matrix as a .txt file (can be thus converted into .csv very easily)
mat2txt, matrix(country_weights) saving(cweights_all.txt) replace


// Generate the residuals/treatment effect data for Germany
// drop the cases_lagged_synth to keep our dataset simple as we create the residuals variable
capture drop cases_lagged_synth
// get the Covid-19 cases predicted by synthetic Germany and save them in a matrix
matrix cases_synth_temp = e(Y_synthetic)

// save predicted synthetic cases in column next to Germany using svmat to convert our matrix to a variable (cases_synth_temp)
replace order = 1 if country_id != `germany_id'
sort order week
svmat cases_synth_temp
// save cases_synth_temp (only synthetic cases for Germany itself) to our final tally of synthetic cases for each synthetic country (the synth_cases variable)
gen synth_cases = cases_synth_temp

// calculate predicted difference (residuals or treatment effect) between Covid-19 cases in synthetic and real Germany
gen treatment_effect = cases_synth_temp - cases_lagged if country_id == `germany_id'
// create a variable containing the residual for Germany (same variable as above but will not be created for each placebo test, this is useful for our later plot overlaying the different residuals when leaving different control countries out)
gen residual_Germany = cases_synth_temp - cases_lagged if country_id == `germany_id'
replace order = 0
// drop cases_synth_temp to clean up (since we already saved the numbers in synth_cases)
drop cases_synth_temp
matrix drop cases_synth_temp


// generate the MSE data
// square the residuals and store them in temporary MSEvar variable
gen MSEvar = treatment_effect^2
// find the mean of the squared residuals for Germany in the time before treatment was implemented to get the pre-treatment MSE
sum MSEvar if country_id == `germany_id' & week < 22
scalar preMSE = r(mean)
// find the mean of the squared residuals for Germany in the time after treatment was implemented to get the post-treatment MSE
sum MSEvar if country_id == `germany_id' & week >= 22
scalar postMSE = r(mean)
// build a matrix with pre-treatment MSE and post-treatment MSE in the two columns (to help us save it later)
matrix MSE1 = (preMSE, postMSE)
// name the columns of the MSE matrix
matrix colnames MSE1 = "pre-MSE" "post-MSE"

// save the rMSPE (root of the mean square predicted error) generated by the synth function in a matrix
matrix rmse = e(RMSPE)
// name the columns of the rmse matrix
matrix colnames rmse = "Tot-RMSE"
// name the row of the rmse matrix (which will be combined with the MSE matrix, naming the row helps us keep a table when we do our permutation-based inference later on)
matrix rownames rmse = "Germany"

// combine the MSE and rmse matrices in a final MSE matrix containing both root mean squared prediction error for the whole period and pre- and post-treatment mean squared error
matrix MSE = rmse, MSE1
// save the MSE matrix as a .txt file (can be thus converted into .csv very easily)
mat2txt, matrix(MSE) saving(MSE_all.txt) replace
// add in the ratio of post-treatment to pre-treatment MSE in the MSErat variable (during permutation-based inference, we will add more ratios to this variable for donor units)
gen MSErat = postMSE/preMSE if country_id == `germany_id'
// drop temporary MSEvar variable to clean up
drop MSEvar

// Here we begin our inference using permutation tests (note that this loop will skip Germany)
// we loop through each country in the EU+ donor pool and generate a synthetic control for it as well, comparing its residuals and MSE ratio to that of the treated unit, Germany
forvalues i = 1(1)30 {

	// skip to next iteration of loop if i represents Germany
	if `i' == `germany_id' {
		continue
	}
	
	// create synthetic donor country from rest of EU+ donor pool + Germany
	synth cases_lagged cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus, trunit(`i') trperiod(22) margin(0.05)

	// generate the residuals/treatment effect data for our current EU+ country
	// drop the cases_lagged_synth to keep our dataset simple as we create the residuals variable	
	capture drop cases_lagged_synth`i'
	matrix cases_synth_temp`i' = e(Y_synthetic)
	
	// save predicted synthetic cases in column next to respective EU+ donor country using svmat to convert our matrix to a variable (cases_synth_temp)
	replace order = 1 if country_id != `i'
	sort order week
	svmat cases_synth_temp`i'
	// save cases_synth_temp`i' (only synthetic cases for current donor country) to our final tally of synthetic cases for each synthetic country (the synth_cases variable)
	replace synth_cases = cases_synth_temp`i' if country_id == `i'

	// calculate predicted difference (residuals or treatment effect) between Covid-19 cases in synthetic and real EU+ donor country
	replace treatment_effect = cases_synth_temp`i' - cases_lagged if country_id == `i'
	replace order = 0
	// drop cases_synth_temp`i' to clean up (since we already saved the numbers in synth_cases)
	drop cases_synth_temp`i'
	matrix drop cases_synth_temp`i'
	
	// generate the MSE data
	// square the residuals and store them in temporary MSEvar variable		
	gen MSEvar = treatment_effect^2
	// find the mean of the squared residuals for EU+ donor country in the time before treatment was implemented to get the pre-treatment MSE
	sum MSEvar if country_id == `i' & week < 22
	scalar preMSE = r(mean)
	// find the mean of the squared residuals for EU+ donor country in the time after treatment was implemented to get the post-treatment MSE
	sum MSEvar if country_id == `i' & week >= 22
	scalar postMSE = r(mean)
	// build a matrix with pre-treatment MSE and post-treatment MSE in the two columns (to help us save it later)
	matrix MSE1 = (preMSE, postMSE)
	// name the columns of the MSE matrix
	matrix colnames MSE1 = "pre-MSE" "post-MSE"

	// save the rMSPE (root of the mean square predicted error) generated by the synth function in a matrix
	matrix rmse = e(RMSPE)
	// name the column of the rmse matrix
	matrix colnames rmse = "Tot-RMSE"
	// find the EU+ donor country name (in arbitrary week) to name the column row
	levelsof country if country_id==`i' & week==3, local(cname) clean
	//scalar cname = country if country_id==`i' & week==3
	// name the row of the rmse matrix (which will be combined with the MSE matrix, naming the row helps us keep a table to compare these residuals with Germany's residuals)
	matrix rownames rmse = `cname'

	// combine the MSE and rmse matrices in a final MSE matrix containing both root mean squared prediction error for the whole period and pre- and post-treatment mean squared error
	matrix MSE = rmse, MSE1
	// append this new donor country MSE matrix to the full MSE .txt file containing in each row pre-, post-, and r-, MSE values for each country in the current analysis (including Germany)
	mat2txt, matrix(MSE) saving(MSE_all.txt) append

	// add in the ratio of post-treatment to pre-treatment MSE in the MSErat variable to compare with other donor country MSE ratios and with the MSE ratio of Germany itself
	replace MSErat = postMSE/preMSE if country_id == `i'
	// drop temporary MSEvar variable to clean up
	drop MSEvar
	
}

// save the resulting dataset as a result dataset
save synthetic_results_full.dta, replace
drop residual_Germany


// plot all residuals from each donor country plus Germany on a graph
sort country week
// seperate treatment effects by country in order to plot them neatly overlaid on a plot
seperate treatment_effect, by(country)

// generate a local macro to contain the line with all EU+ donor residuals as dashed gray lines and Germany's residual as a black solid line
local call
// loop through each country to add dashed gray line to macro
forval j = 1/30 {
	local call `call' line treatment_effect`j' week, lpattern(dash) lcolor(gray) legend(off) || 
}

// add solid black line to macro for Germany's residuals plus configure plot options (adding lines at beginning and end of treatment)
local call `call' line treatment_effect`germany_id' week, ytitle("Residuals of lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) lpattern(solid)  lcolor(black) lwidth(thick) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(dash) lcolor(red)) yline(0, lstyle(foreground) lpattern(solid) lcolor(gs8))

// call twoway function to plot the command contained in the call macro
twoway `call'
// save plot of residuals (as a type of visual inference)
graph save residuals_all, replace
graph export residuals_all.eps, replace


// graph MSE ratios as horizontal bars, again for inference purposes
graph hbar MSErat, over(country, sort(1) label(labsize(vsmall))) ytitle("Ratio of Post- and Pre-9-Euro-Ticket MSPE")
// save plot of MSE ratios
graph save mseratio_all, replace
graph export mseratio_all.eps, replace

// plot residual for Germany by itself
line treatment_effect`germany_id' week, ytitle("Residuals of lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) lpattern(solid)  lcolor(black) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(dash) lcolor(red)) yline(0, lstyle(foreground) lpattern(solid) lcolor(gs8))

graph save synthetic_residual_all, replace
graph export synthetic_residual_all.eps, replace


// plot real Germany vs synthetic Germany
gen synth_cases_germany = cases_lagged11+treatment_effect`germany_id'
twoway line synth_cases_germany week, lpattern(dash) lcolor(gray) || line cases_lagged`germany_id' week, ytitle("Lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) legend(order(1 "Synthetic Germany" 2 "Germany") position(6) rows(1) region(lcolor(black))) lpattern(solid) lcolor(black) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(solid) lcolor(red))

graph save synthetic_germany_all, replace
graph export synthetic_germany_all.eps, replace

drop synth_cases_germany


// save data as synth_results (remember that these do not include Luxembourg)
save synth_results.dta, replace

// clean up and clear
macro drop _all
clear matrix
clear



///////////////////////////////////////////////////////////////////////////////////////////
// ADDITIONAL ROBUSTNESS CHECKS BY MODIFYING THE DONOR GROUP TO EXCLUDE DIFFERENT GROUPS //
///////////////////////////////////////////////////////////////////////////////////////////



// we loop through each country (except Luxembourg which we purposefully excluded from the beginning and Germany which is the treated unit) and leave it out so that it cannot constitute part of synthetic Germany, allowing us to evaluate whether any particular country influences the results too much

// in addition, we perform two special case group exclusions: first, the "Mandates" exclusion which excludes Denmark, Sweden, Finland, and Norway, all the countries that did not have any form of mask mandate, and second, the "Neighbors" exclusion which strengthens our argument for SUTVA by excluding all countries bordering Germany


// loop through each country and the 2 special cases to exclude each country once on it's own and each special case group once as a group
foreach cntry in "Austria" "Belgium" "Bulgaria" "Croatia" "Cyprus" "Czechia" "Denmark" "Estonia" "Finland" "France" "Greece" "Hungary" "Iceland" "Ireland" "Italy" "Latvia" "Liechtenstein" "Lithuania" "Malta" "Netherlands" "Norway" "Poland" "Portugal" "Romania" "Slovakia" "Slovenia" "Spain" "Sweden" "Switzerland" "Mandates" "Neighbors" {
	
	
	display "`cntry'"
	use synth_results.dta
	
	// most of this analysis only leaves one country out but 2 special cases perform our analysis while dropping specific groups of countries		
	if "`cntry'" == "Mandates" {
		// we perform one analysis dropping countries without any type of mask mandate (Denmark, Sweden, Finland, Norway)
		drop if country == "Denmark"
		drop if country == "Sweden"
		drop if country == "Finland"
		drop if country == "Norway"
		// define variable that tells us how many countries we have in the donor pool
		local num_cntry = 26
		display "Excluding Countries without Mask Mandates"
	}
	else if "`cntry'" == "Neighbors" {
		// we perform one analysis dropping countries that neighbour Germany (in order to ensure that no SUTVA violations are occurring)
		drop if country == "Austria"
		drop if country == "Belgium"
		drop if country == "Czechia"
		drop if country == "Denmark"
		drop if country == "France"
		drop if country == "Netherlands"
		drop if country == "Poland"
		// define variable that tells us how many countries we have in the donor pool
		local num_cntry = 23
		display "Excluding Countries without Neighbors"
	}
	else {
		// we perform one analysis where we drop only one country at a time
		drop if country == "`cntry'"
		// define variable that tells us how many countries we have in the donor pool
		local num_cntry = 29
	}

	
	
	// sort countries first by country then by week (such that they get encoded in alphabetical order for aesthetic purposes)
	sort country week

	// set countries to units and weeks to time units
	encode country, generate(country_id2)
	drop country_id
	rename country_id2 country_id
	tsset country_id week

	
	// Find the number by which Germany is coded in the dataset
	gen temp_id = country_id if country == "Germany"
	qui sum temp_id
	local germany_id = `r(max)' 

	// germany_id contains the number coding Germany
	drop temp_id
	

	// create synthetic Germany from EU+ states, with either the individual country or group in the loop having been excluded from that donor pool
	synth cases_lagged cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus, trunit(`germany_id') trperiod(22) margin(0.05) fig
	
	// Uncomment to save graph of synthetic Germany with exclusion, generated directly by the synth function
	//graph save synthetic_germany_no_`cntry', replace
	//graph export synthetic_germany_no_`cntry'.eps, replace
	
	// create a balance table (showing the mean values across our time period for each of the control and target variables for both synthetic and real Germany and the EU+ donor pool)
	
	// loop through each variable used to create synthetic Germany (both controls and targets) and find their mean for entire donor pool
	foreach i in cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus {
		sum `i'
		local mean`i' = r(mean)
	}
	

	// build a matrix containing the average values for each variable for the donor pool
	matrix donor_avg_mat = (`meancases_lagged'\ `meantot_vacc_rate'\ `meanpop_sexratio'\ `meangovstringency'\ `meanpop_density'\ `meanlgdp_cap'\ `meanage_0to19'\ `meanage_20to39'\ `meanage_40to59'\ `meanage_60to79'\ `meanage_80to99'\ `meanage_100plus')
	// name the matrix column
	matrix colnames donor_avg_mat = "Donor Pool Mean"
	// save the balance table generated by the synth function as a matrix containing the variable values for synthetic and real Germany
	matrix balance1 = e(X_balance)
	// combine the matrices containing variable means for synthetic and real Germany and variable means for the donor pool
	matrix balance = balance1, donor_avg_mat
	// save the final balance matrix with all 3 columns as a .txt file (can be thus converted into .csv very easily)
	mat2txt, matrix(balance) saving(balance_no_`cntry'.txt) replace


	// save the weights table generated by the synth function as a matrix containing the weights of each individual donor to synthetic Germany
	matrix country_weights = e(W_weights)
	// save the weights matrix as a .txt file (can be thus converted into .csv very easily)
	mat2txt, matrix(country_weights) saving(cweights_no_`cntry'.txt) replace

	
	// generate the residuals/treatment effect data for Germany
	// drop the cases_lagged_synth to keep our dataset simple as we create the residuals variable
	capture drop cases_lagged_synth
	// get the Covid-19 cases predicted by synthetic Germany and save them in a matrix
	matrix cases_synth_temp = e(Y_synthetic)

	// save predicted synthetic cases in column next to Germany using svmat to convert our matrix to a variable (cases_synth_temp)
	replace order = 1 if country_id != `germany_id'
	sort order week
	svmat cases_synth_temp
	// save cases_synth_temp (only synthetic cases for Germany itself) to our final tally of synthetic cases for each synthetic country (the synth_cases variable)
	// This will be important when we compare MSEs and residuals between countries later on
	replace synth_cases = cases_synth_temp

	// calculate predicted difference (residuals or treatment effect) between Covid-19 cases in synthetic and real Germany
	replace treatment_effect = cases_synth_temp - cases_lagged if country_id == `germany_id'
	
	
	
	// if the country currently being left out belongs to the component countries for our first synthetic Germany, then we save the treatment effect in a variable named residuals which we can use to plot the synthetic Germanys created while leaving each of these countries out overlaid on the same plot
	if "`cntry'" == "Austria" || "`cntry'" == "Finland" || "`cntry'" == "France" || "`cntry'" == "Italy" || "`cntry'" == "Malta" {
		gen residual_no_`cntry' = cases_synth_temp - cases_lagged if country_id == `germany_id'
		// generate a variable that we can use to merge into the synthetic_results_full dataset later
		gen residual_name = country
		replace residual_name = "`cntry'" if country_id == `germany_id'
		save synth_results_no_`cntry'.dta, replace
		use synthetic_results_full.dta
		// create corresponding merging variable (just a copy of the country variable)
		gen residual_name = country
		// merge
		merge 1:1 residual_name week using synth_results_no_`cntry'.dta, keepusing(residual_no_`cntry') nogenerate
		// drop the merging variable so we can regenerate it later if we need it again
		drop residual_name
		save synthetic_results_full.dta, replace
		use synth_results_no_`cntry'.dta
		
	}
	replace order = 0
	// drop cases_synth_temp to clean up (since we already saved the numbers in synth_cases)
	matrix drop cases_synth_temp
	drop cases_synth_temp


	// generate the MSE data
	// square the residuals and store them in temporary MSEvar variable
	gen MSEvar = treatment_effect^2
	// find the mean of the squared residuals for Germany in the time before treatment was implemented to get the pre-treatment MSE
	sum MSEvar if country_id == `germany_id' & week < 22
	scalar preMSE = r(mean)
	// find the mean of the squared residuals for Germany in the time after treatment was implemented to get the post-treatment MSE
	sum MSEvar if country_id == `germany_id' & week >= 22
	scalar postMSE = r(mean)
	// build a matrix with pre-treatment MSE and post-treatment MSE in the two columns (to help us save it later)
	matrix MSE1 = (preMSE, postMSE)
	// name the columns of the MSE matrix
	matrix colnames MSE1 = "pre-MSE" "post-MSE"

	// save the rMSPE (root of the mean square predicted error) generated by the synth function in a matrix
	matrix rmse = e(RMSPE)
	// name the columns of the rmse matrix
	matrix colnames rmse = "Tot-RMSE"
	// name the row of the rmse matrix (which will be combined with the MSE matrix, naming the row helps us keep a table when we do our permutation-based inference later on)
	matrix rownames rmse = "Germany"

	// combine the MSE and rmse matrices in a final MSE matrix containing both root mean squared prediction error for the whole period and pre- and post-treatment mean squared error
	matrix MSE = rmse, MSE1
	// save the MSE matrix as a .txt file (can be thus converted into .csv very easily)
	mat2txt, matrix(MSE) saving(MSE_no_`cntry'.txt) replace
	
	// add in the ratio of post-treatment to pre-treatment MSE in the MSErat variable (during permutation-based inference, we will add more ratios to this variable for donor units)
	replace MSErat = postMSE/preMSE if country_id == `germany_id'
	// drop temporary MSEvar variable to clean up
	drop MSEvar



	// Here we begin our inference using permutation tests (note that this list contains all but the left out country)
	// we loop through each country in our list and generate a synthetic control for it as well, finally comparing it's residuals and MSE ratio to that of the treated unit, Germany
	forvalues i = 1(1)`num_cntry' {
	
		// skip to next iteration of loop if i represents Germany
		if `i' == `germany_id' {
			continue
		}


		// create synthetic donor country from rest of donor pool + germany
		synth cases_lagged cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus, trunit(`i') trperiod(22) margin(0.05)


		// generate the residuals/treatment effect data for our current donor country
		// drop the cases_lagged_synth to keep our dataset simple as we create the residuals variable
		capture drop cases_lagged_synth`i'
		matrix cases_synth_temp`i' = e(Y_synthetic)
	
		// save predicted synthetic cases in column next to respective donor country using svmat to convert our matrix to a variable (cases_synth_temp)
		replace order = 1 if country_id != `i'
		sort order week
		svmat cases_synth_temp`i'
		// save cases_synth_temp`i' (only synthetic cases for current donor country) to our final tally of synthetic cases for each synthetic country (the synth_cases variable)
		replace synth_cases = cases_synth_temp`i' if country_id == `i'

		// calculate predicted difference (residuals or treatment effect) between Covid-19 cases in synthetic and real donor country
		replace treatment_effect = cases_synth_temp`i' - cases_lagged if country_id == `i'
		replace order = 0
		// drop cases_synth_temp`i' to clean up (since we already saved the numbers in synth_cases)
		matrix drop cases_synth_temp`i'
		drop cases_synth_temp`i'


		// generate the MSE data
		// square the residuals and store them in temporary MSEvar variable		
		gen MSEvar = treatment_effect^2
		// find the mean of the squared residuals for donor country in the time before treatment was implemented to get the pre-treatment MSE
		sum MSEvar if country_id == `i' & week < 22
		scalar preMSE = r(mean)
		// find the mean of the squared residuals for donor country in the time after treatment was implemented to get the post-treatment MSE
		sum MSEvar if country_id == `i' & week >= 22
		scalar postMSE = r(mean)
		// build a matrix with pre-treatment MSE and post-treatment MSE in the two columns (to help us save it later)
		matrix MSE1 = (preMSE, postMSE)
		// name the columns of the MSE matrix
		matrix colnames MSE1 = "pre-MSE" "post-MSE"

		// save the rMSPE (root of the mean square predicted error) generated by the synth function in a matrix
		matrix rmse = e(RMSPE)
		// name the columns of the rmse matrix
		matrix colnames rmse = "Tot-RMSE"
		// find the donor country name (in arbitrary week) to name the column row
		levelsof country if country_id==`i' & week==3, local(cname) clean
		// name the row of the rmse matrix (which will be combined with the MSE matrix, naming the row helps us keep a table to compare the donor residuals with Germany's residuals)
		matrix rownames rmse = `cname'

		// combine the MSE and rmse matrices in a final MSE matrix containing both root mean squared prediction error for the whole period and pre- and post-treatment mean squared error
		matrix MSE = rmse, MSE1
		// append this new donor country MSE matrix to the full MSE .txt file containing in each row pre-, post-, and r-, MSE values for each country (including Germany)
		mat2txt, matrix(MSE) saving(MSE_no_`cntry'.txt) append
	
		// add in the ratio of post-treatment to pre-treatment MSE in the MSErat variable (allowing us to compare the donor country MSE ratio to Germany's MSE ratio)
		replace MSErat = postMSE/preMSE if country_id == `i'
		// drop temporary MSEvar variable to clean up
		drop MSEvar
	}

	// plot all residuals from each donor country plus Germany on a graph

	// loop through each country and drop their individual treatment effect variables (for cleanup and since we want to seperate them more effectively)
	sort country week
	forval j = 1/`num_cntry' {
		capture drop treatment_effect`j'
	}
	seperate treatment_effect, by(country)
	
	// generate a local macro to contain the line with all donor residuals as dashed gray lines and Germany's residual as a black solid line
	local call
	// loop through each country to add dashed gray line to macro
	forval j = 1/`num_cntry' {
		local call `call' line treatment_effect`j' week, lpattern(dash) lcolor(gray) legend(off) || 
	}

	// add solid black line to macro for original synthetic Germany's residuals (excluding none of the donor pool) plus configure plot options (adding lines at beginning and end of treatment)
	local call `call' line treatment_effect`germany_id' week, ytitle("Residuals of lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) lpattern(solid)  lcolor(black) lwidth(thick) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(solid) lcolor(red)) yline(0, lstyle(foreground) lpattern(solid) lcolor(gs8))

	


// call twoway function to plot the command contained in the call macro
	twoway `call'
	// save plot of residuals (as a type of visual inference)
	graph save residuals_no_`cntry', replace
	graph export residuals_no_`cntry'.eps, replace

	// graph MSE ratios as horizontal bars, again for inference purposes
	graph hbar MSErat, over(country, sort(1) label(labsize(vsmall))) ytitle("Ratio of Post- and Pre-9-Euro-Ticket MSPE")
	// save plot of MSE ratios
	graph save mseratio_no_`cntry', replace
	graph export mseratio_no_`cntry'.eps, replace
	
	
	gen synth_cases_germany = cases_lagged11+treatment_effect`germany_id'
	// plot real Germany vs synthetic Germany without 'cntry'
	twoway line synth_cases_germany week, lpattern(dash) lcolor(gray) || line cases_lagged11 week, ytitle("Lagged COVID-19 Cases per Capita") xtitle(		"Calendar Week") xlabel(0(5)40) legend(order(1 "Synthetic Germany" 2 "Germany") position(6) rows(1) region(lcolor(black))) lpattern(solid) lcolor(black) 		xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(solid) lcolor(red))

	graph save synthetic_germany_no_`cntry', replace
	graph export synthetic_germany_no_`cntry'.eps, replace
	
	drop synth_cases_germany

	// save the modified dataset for this particular configuration of donor group
	save synth_results_no_`cntry'.dta, replace
	
	// clean up and clear
	clear matrix
	clear

}


// now we use the synthetic_results_full which contains the residuals for the original synthetic versus actual Germany analysis to plot 
use synthetic_results_full.dta

// generate a local macro to contain the line with all synthetic Germany residuals (excluding the original component countries one at a time) as dashed gray lines and Germany's residual as a black solid line
local call
// loop through each country to add dashed gray line to macro
foreach cntry in "Austria" "Finland" "France" "Italy" "Malta" {
	local call `call' line residual_no_`cntry' week, lpattern(solid) lcolor(gray) legend(off) || 
}

// add solid black line to macro for Germany's original residuals plus configure plot options (adding lines at beginning and end of treatment)
local call `call' line residual_Germany week, ytitle("Residuals of lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) lpattern(dash)  lcolor(black) lwidth(thick) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(solid) lcolor(red)) yline(0, lstyle(foreground) lpattern(solid) lcolor(gs8))

// call twoway function to plot the command contained in the call macro
twoway `call'
// save plot of overlaid synthetic Germany residuals (as main plot for leave-one-out analysis)
graph save residuals_overlaid_loo, replace
graph export residuals_overlaid_loo.eps, replace



clear matrix
macro drop _all
clear


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ADDITIONAL TEST USING MALTA'S FREE PUBLIC TRANSPORTATION PROGRAM ON OCTOBER 1st 2022 AS ANOTHER HYPOTHETICAL TREATMENT //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// reset file location to grab correct files and start over anew

cd ..
cd "Compiled Files"
use global_covid.dta
cd ..
cd "Result Files"


// plot cases by state overlaid
// Note that week 39 was when the program was started
seperate cases_lagged, by(country)
gen order = 0


// drop Luxembourg since Luxembourg also had similar public transport subsidy
// drop Germany since we are looking at Malta and Germany had public transport subsidy before Malta
drop if country == "Luxembourg"
drop if country == "Germany"


// set countries to units and weeks to time units for panel data use
encode country, generate(country_id)
tsset country_id week


// Find the number by which Malta is coded in the dataset
gen temp_id = country_id if country == "Malta"
qui sum temp_id
local malta_id = `r(max)' 

// malta_id contains the number coding Malta
drop temp_id


// create synthetic Malta from EU+ donor countries, week 39 is treatment begin for Malta (i.e., week of October 1st 2022)
synth cases_lagged cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus, trunit(`malta_id') trperiod(39) margin(0.05) fig

// Uncomment to save graph of synthetic versus actual lagged Covid-19 cases in Malta, directly generated by synth function
//graph save synthetic_malta, replace
//graph export synthetic_malta.eps, replace

// generate balance table (follows same procedure as above)
foreach i in cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus {
	sum `i'
	local mean`i' = r(mean)
}

// create a matrix containing average numbers for each of the control and target variables
matrix donor_avg_mat = (`meancases_lagged'\ `meantot_vacc_rate'\ `meanpop_sexratio'\ `meangovstringency'\ `meanpop_density'\ `meanlgdp_cap'\ `meanage_0to19'\ `meanage_20to39'\ `meanage_40to59'\ `meanage_60to79'\ `meanage_80to99'\ `meanage_100plus')
matrix colnames donor_avg_mat = "Donor Pool Mean"
// extract and save in new matrix the mean values for the control and target variables for synthetic and real Malta
matrix balance1 = e(X_balance)
matrix balance = balance1, donor_avg_mat
// save balance matrix as .txt table (easily convertible and usable as de facto .csv file)
mat2txt, matrix(balance) saving(balance_malta.txt) replace

// save matrix containing the relative weights of constituent donor countries for synthetic Malta
matrix country_weights = e(W_weights)
mat2txt, matrix(country_weights) saving(cweights_malta.txt) replace


// get the predicted synthetic Covid-19 cases from synthetic Malta
capture drop cases_lagged_synth
matrix cases_synth_temp = e(Y_synthetic)

// save predicted synthetic cases in column next to Malta using svmat function
replace order = 1 if country_id != `malta_id'
sort order week
svmat cases_synth_temp

// calculate predicted difference (residuals)
gen synth_cases = cases_synth_temp
gen treatment_effect = cases_synth_temp - cases_lagged if country_id == `malta_id'
replace order = 0
drop cases_synth_temp
matrix drop cases_synth_temp

// generate MSE data as above
gen MSEvar = treatment_effect^2
sum MSEvar if country_id == `malta_id' & week < 39
scalar preMSE = r(mean)
sum MSEvar if country_id == `malta_id' & week >= 39
scalar postMSE = r(mean)

matrix MSE1 = (preMSE, postMSE)
matrix colnames MSE1 = "pre-MSE" "post-MSE"

matrix rmse = e(RMSPE)
matrix colnames rmse = "Tot-RMSE"
matrix rownames rmse = "Malta"

matrix MSE = rmse, MSE1
mat2txt, matrix(MSE) saving(MSE_malta.txt) replace

// generate MSE ratio for Malta (to be compared with other placebo tests on synthetic donor countries)
gen MSErat = postMSE/preMSE if country_id == `malta_id'
drop MSEvar

// inference using permutation tests (note that this list does not include Malta itself)
forvalues i = 1(1)29 {
	
	if `i' == `malta_id' {
		continue
	}
	
	// generate synthetic donor country composed of other donor pool countries + Malta
	synth cases_lagged cases_lagged tot_vacc_rate pop_sexratio govstringency pop_density lgdp_cap age_0to19 age_20to39 age_40to59 age_60to79 age_80to99 age_100plus, trunit(`i') trperiod(39) margin(0.05)

	// get the predicted synthetic covid-19 cases for current donor country
	capture drop cases_lagged_synth`i'
	matrix cases_synth_temp`i' = e(Y_synthetic)
	
	// save predicted synthetic cases in column next to the respective donor country
	replace order = 1 if country_id != `i'
	sort order week
	svmat cases_synth_temp`i'


	// calculate predicted difference (residuals) for donor country
	replace synth_cases = cases_synth_temp`i' if country_id == `i'
	replace treatment_effect = cases_synth_temp`i' - cases_lagged if country_id == `i'
	replace order = 0
	drop cases_synth_temp`i'
	matrix drop cases_synth_temp`i'
	
	// generate MSE data for donor country
	gen MSEvar = treatment_effect^2
	sum MSEvar if country_id == `i' & week < 39
	scalar preMSE = r(mean)
	sum MSEvar if country_id == `i' & week >= 39
	scalar postMSE = r(mean)
	matrix MSE1 = (preMSE, postMSE)
	
	matrix colnames MSE1 = "pre-MSE" "post-MSE"
	matrix rmse = e(RMSPE)
	matrix colnames rmse = "Tot-RMSE"
	
	levelsof country if country_id==`i' & week==3, local(cname) clean
	matrix rownames rmse = `cname'
	
	matrix MSE = rmse, MSE1
	mat2txt, matrix(MSE) saving(MSE_malta.txt) append
	
	// add on to MSE ratio variable with this donor country's MSE ratio to be compared to other donor country placebo runs and Malta's MSE ratio as well
	replace MSErat = postMSE/preMSE if country_id == `i'
	drop MSEvar

	
}

// plot placebo treatment effects for all donor countries and actual treatment effect for Malta (for visual inference purposes)
sort country week
seperate treatment_effect, by(country)
seperate MSErat, by(country)


local call
forval j = 1/29 {
	local call `call' line treatment_effect`j' week, lpattern(dash) lcolor(gray) legend(off) || 
}

local call `call' line treatment_effect`malta_id' week, ytitle("Residuals of lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) lpattern(solid)  lcolor(black) lwidth(thick) xline(39, lstyle(foreground) lpattern(solid) lcolor(green)) yline(0, lstyle(foreground) lpattern(solid) lcolor(gs8))
twoway `call'
graph save residuals_malta, replace
graph export residuals_malta.eps, replace

// graph MSE ratios
graph hbar MSErat, over(country, sort(1) label(labsize(vsmall))) ytitle("Ratio of Post- and Pre-October MSPE") 

graph save mseratio_malta, replace
graph export mseratio_malta.eps, replace


// plot real Malta vs synthetic Malta
gen synth_cases_malta = cases_lagged`malta_id'+treatment_effect`malta_id'
twoway line synth_cases_malta week, lpattern(dash) lcolor(gray) || line cases_lagged`malta_id' week, ytitle("Lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) legend(order(1 "Synthetic Malta" 2 "Malta") position(6) rows(1) region(lcolor(black))) lpattern(solid) lcolor(black) xline(22, lstyle(foreground) lpattern(solid) lcolor(green)) xline(34, lstyle(foreground) lpattern(solid) lcolor(red))

graph save synthetic_malta, replace
graph export synthetic_malta.eps, replace

drop synth_cases_malta


// plot residual for Malta by itself
line treatment_effect`malta_id' week, ytitle("Residuals of lagged COVID-19 Cases per Capita") xtitle("Calendar Week") xlabel(0(5)40) lpattern(solid) lcolor(black) xline(39, lstyle(foreground) lpattern(solid) lcolor(green)) yline(0, lstyle(foreground) lpattern(solid) lcolor(gs8))

graph save synthetic_residual_malta, replace
graph export synthetic_residual_malta.eps, replace

// save data as synth_results (remember that these do not include Luxembourg or Germany)
save synth_results_malta.dta, replace

cd ..


////////////////
//END OF FILE //
////////////////


