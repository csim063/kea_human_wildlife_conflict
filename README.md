# Kennedy, Perry & Simpkins (2020). Additive influence of human-wildlife conflict and introduced mammalian predation on the population dynamics of kea (*Nestor notabilis*)
Repository contains the data and code needed to redo the analyses described in the manuscript, as well as an RMarkdown file used to recreate the manuscript (and figures) itself.

# Abstract
Human populations are continuing to expand into previously wild areas, bringing humans and wildlife into conflict. Human-wildlife conflicts often directly result in wildlife mortality events, termed human-induced mortality (HM). While HM is a well known phenomena, its overall impacts at the population-level has not been well studied. HM may interact with other sources of mortality, such as invasive predators, exacerbating both beyond the levels accounted for in conservation management plans. Our study aimed to explore the population-level impacts of changes in HM intensity, and how these impacts would interact with a major known source of mortality. To achieve these aims we developed a stage-based population model to simulate the population dynamics of kea (*Nestor notabilis*), an endangered parrot endemic to New Zealand. We used this model to run multiple scenarios with differing intensity levels of invasive mammalian predation and human-induced mortality. Mammalian predation had the most pronounced impact on kea populations. With unmanaged predation resulting in rapid extinction of the population. HM had a far smaller impact, with the current rate of HM not severely affecting kea population dynamics. However, when HM grew continually over time, simulating increased human populations, the kea population showed significant decreases in population size and extinction risk over time, and this was exacerbated by mammalian predation, even at currently managed levels. These results clearly show that while HM may not be an immediately pressing threat, if left unmanaged it can rapidly become a major issue to conservation management.

# Repository design
## Data
Most of the data used in the study were values sourced from literature or Department of Conservation databases, and as such are provided as values in the literature itself and are not stored separately in this repository. As the time taken to run the model can be fairly substantial, however, a number of result .RData files are stored within `Output` to allow for a more rapid building of the RMarkdown manuscript. The `External_data` folder contains data relating to kea siting locations and a spatial polygon of New Zealand's coastlines used to develop Figure 1 in the manuscript.

## Scripts
The `src` folder contains all the scripts used to create, run, and analyse the population viability model.
- `Analyse_results.R`: Analyses the data output from the PVA model run using `Run_analysis.R`
- `Input_variables.R`: Used to set all the initial and constant values used in the PVA model
- `PVA_Function.R`: Creates the function used in the ode() function. This function is written with all the variables matching those required by the deSolve functions. It also returns a list of the individuals in each age class, as well as the total population. **The main script which forms the model**
- `Run_analysis_function.R`: Creates a function to run the PVA model
- `Run_sensitivity_analysis.R`: Used to run the PVA model for the sensitivity analysis

The model is run using the `Run_analysis.R` found in the root directory.

## Manuscript
The `manuscript` folder contains `kea_PVA_manuscript.Rmd` and its outputs (.tex and .pdf), which was used to create the submitted manuscript. The `figures` folder contains all the figures not produced by the .Rmd file itself. `kea_PVA_manuscript_files/figure-latex` contains all the figures created by the .Rmd file. `References` contains `master.bib` which is the Bibtex library used to generate all the citations in the manuscript. `manuscript_style_template.tex` is a file used to style the manuscript .pdf output.

# Further information
If you require any additional information or have any issues please describe them in the Issues tab of this repository or email simpkinscraig063@gmail.com