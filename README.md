# public-transit-CoV Repository

## Description
This is the repository containing all the necessary code and data to reproduce the 2024 paper "Ticket to Ride the Covid Train? The Effect of Public Transportation Subsidies on COVID-19". It also contains the results of the analysis, intermediate .dta files used for compiling these results, and the final figures (slightly modified post-generation).

### Data
Most data required for the analysis are directly sourced from the respective public github repositories. COVID-19 case numbers are directly taken from the global confirmed time series data provided by COVID-19 Data Repository by the Center for Systems Science and Engineering at Johns Hopkins University (Dong et al., 2020). The Oxford COVID-19 Government Response Tracker provides a daily index constructed to reflect the severity of COVID-19 containment policies (Hale et al., 2021). Vaccination rates are sourced from the "Our World in Data" COVID-19 dataset (Mathieu et al., 2021).

I was unable to directly link to online datasources for control data from the World Population Prospects survey (providing demographic information) nor from the National Accounts Main Aggregates Database (providing national GDP) (United Nations, 2022; United Nations Department of Economic and Social Affairs, 2022). Corresponding csv files were downloaded instead and placed inside the "CSV Data" folder. Given that this is a public repository, the execution code file directly links to the raw data in this repo, making it unneccessary to download this data on your own.

### Compiled Files and Result Files
Compiled Files contains 2 .dta datasets that represent the data after processing and before performing any systematic analysis. Result Files contains tables, .dta files, and figures (in both .gph and .eps format) that result from running the .do file.

### Code
The annotated_full.do file contains the full analysis seen in the paper, from data processing to synthetic control. Please note our use of the synth package first published by Abadie, Diamond, & Hainmueller (2010). Further details regarding this package can be found at the following url: https://web.stanford.edu/~jhain/synthpage.html.

## Getting Started

### Prerequisites
Stata 18 was used to execute this code. A internet connection is required while running the .do file. The synth, outreg2, and mat2txt packages are installed during execution (can be commented out in the .do in lines 449, 483, and 484).

### Installation
Download the annotated_full.do file. Create a folder called "Compiled Files" and a folder called "Result Files" within the same folder as the .do file. Alternatively, you can just clone this repository.

### Execution
Execute the annoated_full.do file in Stata. Intermediate data tables will be saved in the "Compiled Files" folder, resulting tables and figures will be saved in the "Result Files" folder.

## Authors
Christian Metzner
@cjm972

## License
This project is licensed under the GNU General Public License v3.0. For more details, refer to the LICENSE.md file.

## Citations

Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic Control Methods for comparative case studies: Estimating the effect of California’s tobacco control program. Journal of the American Statistical Association, 105(490), 493–505.

Dong, E., Du, H., & Gardner, L. (2020). An interactive web-based dashboard to track COVID-19 in real time. The Lancet Infectious Diseases, 20(5), 533–534. https://doi.org/10.1016/s1473-3099(20)30120-1

Hale, T., Angrist, N., Goldszmidt, R., Kira, B., Petherick, A., Phillips, T., Webster, S., Cameron-Blake, E., Hallas, L., Majumdar, S., & Tatlow, H. (2021). A global panel database of pandemic policies (Oxford Covid-19 Government response tracker). Nature Human Behaviour, 5(4), 529–538. https://doi.org/10.1038/s41562-021-01079-8

Mathieu, E., Ritchie, H., Ortiz-Ospina, E., Roser, M., Hasell, J., Appel, C., Giattino, C., & Rodés-Guirao, L. (2021). A global database of COVID-19 vaccinations. Nature Human Behaviour, 5(7), 947–953. https://doi.org/10.1038/s41562-021-01122-8

United Nations. (2022). World Population Prospects. United Nations. Retrieved November 13, 2022, from https://population.un.org/wpp/Download/Standard/CSV/ 

United Nations Department of Economic and Social Affairs. (2022). Basic data selection. United Nations. https://unstats.un.org/unsd/snaama/Basic




