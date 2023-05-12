## "Predicting Covid-19 infections using multi-layer centrality measures"
Author: **Christine Hedde-von Westernhagen, May 2023**

### Thesis project for the MSBBSS Master's program 
*Department of Methodology and Statistics, Utrecht University, 2022-2023. Supervisors: Javier Garcia-Bernardo and Ayoub Bagheri.*

##### Abstract

<img src="https://raw.githubusercontent.com/christine-hvw/thesis_disclosed/example_analyses/analyses/3dplot_dummy.png" align="right" width="350px">

*One of the most pressing problems of the recent past has been to understand the spread of the Covid-19 virus within and across populations. A popular strategy in investigating spreading phenomena have been network models, as they are able to represent complex social interactions. While previous studies have shown that network centrality measures are able to identify influential spreaders, it is not yet known if they can also be used to predict individual infection risks. However, information about individual risks is vital for the design of targeted interventions, and can give citizens more autonomy over their decisions. This study therefore focuses on the predictive abilities of centrality measures for individual infections. I investigate this issue leveraging the framework of multi-layer networks, which allows to explicitly consider human interaction to take place in different contexts simultaneously. Drawing on large-scale administrative data from the Netherlands, I find that existing centrality measures offer good predictions of relative infection risks, but are only weakly informative about the timing of individual infections. I therefore introduce a new Degree-based multi-layer centrality measure that takes into account the transmission rate of Covid-19 in different contexts. The new measure shows similar predictive performance as existing centrality measures, indicating that centrality measures alone are limited in their potential to predict the timing of infections, but could be used to complement other prediction approaches.*

---

***Please note:***

- This repository is a publicly accessible version of my thesis project. It serves as a documentation of the performed analyses, so they can be reproduced if necessary. It also stores the final version of the manuscript. There is NO DATA stored in this repository due to privacy regulations.

- To reproduce results shown in the thesis manuscript (`thesis_chw_v080523.pdf`), all code in the file `analyses.Rmd` has to be executed after obtaining access to the original data ([info below](#data-access)).

- The code, data and output on the branch [`example_analyses`](https://github.com/christine-hvw/thesis_disclosed/tree/example_analyses) served experimentation purposes in earlier stages of the project.


---


### Repository Contents

| Content                 | Description                                                                     |
| ---------------------   | ---------------------------------------------------------------------------     |
| 📂`data_processed/`     | stores intermediary data objects (*empty*)                                      |
| 📂`results/`            | stores results of analyses as displayed in `thesis_chw_v080523.pdf` (*empty*)   |
| `analyses.Rmd`          | script describing and executing all data preparation and analyses               |
| `aux_functions.R`       | script with auxiliary functions used in `analyses.Rmd`                          |
| `session_info.txt`      | file containing R session info including package versions                       |
| `thesis_chw_v080523.pdf`| manuscript of thesis as submitted on 8 May 2023                                 |


<a id="data-access"></a>

### Data Access

*Disclaimer*: Results displayed in `thesis_chw_v080523.pdf` are based on calculations by the author using non-public micro-data from Statistics Netherlands (CBS). Under certain conditions, these micro-data are accessible for statistical and scientific research. For further information contact microdata@cbs.nl. Permission to use the data for the presented analyses has been granted by the Ethical Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University on November 8, 2022 (case numbers 22-1886, 22-1887, 22-1888).

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
