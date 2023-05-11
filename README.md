## "Predicting Covid-19 infections using multi-layer centrality measures"
**Christine Hedde-von Westernhagen, May 2023**

### Thesis project for the MSBBSS Master's program 
*Department of Methodology and Statistics, Utrecht University, 2022-2023. Supervisors: Javier Garcia-Bernardo and Ayoub Bagheri.*

---

This repository is a publicly accessible version of my thesis project. It serves as a documentation of the performed analyses, so they can be reproduced if necessary. It also stores the final version of the manuscript. There is NO DATA stored in this repository due to privacy regulations.

To reproduce results shown in the thesis manuscript (`thesis_chw_v080523.pdf`), all code in the file `analyses.Rmd` has to be executed after obtaining access to the original data ([info below](#data-access)).

The code, data and output on the branch [`example_analyses`](link-to-branch) served experimentation purposes in earlier stages of the project.

![grafik](https://github.com/christine-hvw/thesis_disclosed/blob/example-analyses/analyses/3dplot_dummy.png?raw=true)

---

### Repository Contents

| Content                 | Description                                                                     |
| ---------------------   | ---------------------------------------------------------------------------     |
| 📂`data_processed/`     | stores intermediary data objects (*empty*)                                      |
| 📂`results/`            | stores results of analyses  as displayed in `thesis_chw_v080523.pdf` (*empty*)  |
| `analyses.Rmd`          | script describing and executing all data preparation and analyses               |
| `aux_functions.R`       | script with auxiliary functions used in `analyses.Rmd`                          |
| `thesis_chw_v080523.pdf`| manuscript of thesis as submitted on 8 May 2023                                 |


### Data Access {#data-access}

*Disclaimer*: Results are based on calculations by the author using non-public micro-data from Statistics Netherlands (CBS). Under certain conditions, these micro-data are accessible for statistical and scientific research. For further information contact microdata@cbs.nl.

**CBS datasets used**:

Network data ([documentation](https://www.cbs.nl/nl-nl/onze-diensten/maatwerk-en-microdata/microdata-zelf-onderzoek-doen/microdatabestanden/pn-a-persons-netwerk-in-the-netherlands)):

- PersNw2018_v1.0_links_familie.csv
- PersNw2018_v1.0_links_huishouden.csv
- PersNw2018_v1.0_links_werk.csv
- PersNw2018_v1.0_linktype_labels.csv

Schools and their locations ([documentation](https://www.cbs.nl/nl-nl/onze-diensten/maatwerk-en-microdata/microdata-zelf-onderzoek-doen/microdatabestanden/brinadressen-locatie-van-onderwijsinstellingen)):

- BRINADRESSEN2020V1.sav

Registered students ([documentation](https://www.cbs.nl/nl-nl/onze-diensten/maatwerk-en-microdata/microdata-zelf-onderzoek-doen/microdatabestanden/inschrwpotab-inschrijvingen-in-het-basisonderwijs)):

- INSCHRWPOTAB2020V2.sav

Covid tests: 

- CoronIT_GGD_testdata_20210921.sav

---

### Contact

For questions or comments please contact me via c.hedde-vonwesternahagen@students.uu.nl, or open an [issue](https://github.com/christine-hvw/thesis_disclosed/issues) in this repository.
