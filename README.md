# Intorduction
Codes and data of S. Shibasaki, M. Mobilia and S. Mitri (2020) [link](https://arxiv.org/abs/2007.12090).
Here, you can access the programming codes (in Python and C language) and csv files used in the above manuscript. However, we cannot depositi row data (i.e., dynamics in each simulation) due to the limitation of the file sizes. Please contaact Shota Shibasaki if you need such row data.

To run C language codes, you need to have "Mt.h" file.


# mono_extinction
This folder contains a python code to plot Fig 2 in the main text (Fig2_drawer.py) and data files in csv used in Fig2. The mono-culture simulation is implemented using "main_2sp_simulation.c" in folder "interaction_analysis/main_text_scenario1", but one needs small modification such that initial species 2's abundance is 0. Three folders (death2 , death 5, death10) contain csv files that show the abundances of species 1, resource, and toxin as well as the enviornmental condition at the end of each simulation at varisou switching rate. These files are used to show distibutions of species 1's abundance at the end of simulations.

# interaction_analysis
This folder contains two folders: "main_text_scenareio1" and "appendix". Each folder contains codes for simulations and plotting results, as well as summary data used in the main text or appendix.

##  main_text
The main text analyzes the environmental switching scenario 1: switching only resource supplies. In this folder, "main_2sp_simulation.c" runs two-species dynamics where species 1 grows faster than species 2 but the rest parameter values are identical. Then, you will pobtain following folders.

To plots the Figs. 3 and 4 in the main text, use "main_figures.py".

Each folder named "deathDelta", where Delta represents toxin sensitivity, has two folders: "OneConsumer" and "TwoConsumer". In "OneConsumer", you can see the probabilities that species 1 goes extinct in mono-culture over the switching rate. In "TwoConsumer", you will obtain the probabilities that both species 1 and 2 go extinct ("TwoConsumer_BothExtinction_model1_competitor1.csv"), that two species coexist ("TwoConsumer_Coexistence_model1_competitor1.csv"), that species 2 excludes species 1 ("TwoConsumer_CompetitiveExclusion_model1_competitor1.csv"), and that species 1 goes extinct in co-culture ("TwoConsumer_Extinction_model1_competitor1.csv"). In each file each row corresponds to the different species 2 (ratio of maximum growth rates mu1/mu2=1.1, 1.2, ..., 1.5) while each column represents the different switching rate. Note that we used only the first row in the manuscript.

In the folder of "Correlation", you can access the summury of difference in extinctionprobabilities of species 1 in the three different switching scenarios and four different resource supplies in the switching scenario 1.

## appendix
In this folder, you can see codes and data used in Appendices 1, 2, 3, 4, and 5. 
By using the codes in deiveristy_analysis, you can plot figures in Appendix 6.
To plot the figures in appoendix, use "appendix_figure.py" or "main_figures.py".

### Appendix3_constant_env
In this folder, you see the results when there is no envio0rnmental switching. For scenario 1, we show the cases when the resource supply is fixed as scarce, mean or abundant. For scenarios 2 and 3, we show only the cases of mild environments (scarce toxin supply for sceanrio 2, and abundant resource and scarce toxin supllies in scenario 3).

In each scenario, you can see the probabilities (i) that species 1 goes extinc in mono-culture ("OneConsumer_Extinction_), (ii) that species 1 goes extinct in co-culture with species 2 ("TwoConsumer_ExtinctionofSpeciess1"), (iii) of diffrence from (i) to (ii) ("Diff_extinction"), (iv) that both species 1 and 2 go extinct ("BothExtinction"), and (v) that species 2 excludes species 1 (ComeptitiveExclusion). 

### Appendix4_changing_scenario
In this scenario, you can see the results of different environmental switching scenarios: scenario2 (switching toxin supplies) and scenario3 (switching both resource and toxin supplies). In each scenario, you can see (i) extinction probabilities of species 1 in mono-culture, (ii) probabilities that both species go extinct, (iii) probabilities that two species xoexist, (iv) probabilities that species 2 excludes species 1, and (v)species 1's extinction probabiliries in co-culture over the switching rate.

Running simulations and plotting results in these scenarios can be done by modifying "main_2sp_simulation.c" and "main_figures.py".

### Appendix5_change_supply
In this folder, you can see the results of different resource supplies in scernario 1: more abundant resource supply (large_xmax), less abundant resource supply (small_xmax), more scarce resource supply (small_xmin), and less scarce resource supply (large_xmin). In each scenario, you can see (i) extinction probabilities of species 1 in mono-culture, (ii) probabilities that both species go extinct, (iii) probabilities that two species xoexist, (iv) probabilities that species 2 excludes species 1, and (v)species 1's extinction probabiliries in co-culture over the switching rate.

## OtherFlucutuations
This folder contains two results when we implement environment fluctuations other than symmetric switching between two environmental states (harsh and mild): asymmetric switching and cyclic flucutuaion. In "Asymmetric" folder, we can see two scenarios: the harsh environments tend to continues longer (see folder 1) and the mild environment continues longer (folder 2). The containts of each folder are similar to Appendix4. In "FourState" folder we can see the results when we have four environmental states and the enviornmental condition flucutuates such that 1->2->3->4->1->... with a rate \nu.

Running simulations and plotting results in these scenarios can be done by modifying "main_2sp_simulation.c" and "main_figures.py".



# diversity_analysis
In this file, you can see four programming codes (incl. "MT.h") and two data folders (beta_diversity and Species_richness). Each folder contains the python codes that plots the results 

"parm_geenrators.c" generates parameter values of each species in each run. You need to run this code before running the below codes.

"main_diversity.c" performs Glispies algorithm with given numbers of speices and compounds in a system. You will obtain the states of system (amounts of each compopund and population sizes) at the end of each run. In default, this codes run 100  communities with 100 replicate given the number of species, mean toxin sensitivities, and the switching rate.

"competitive_prob.py" analyzes the probability of competitive exclusion (i.e., probability that inferior species in the mean feild excludes the superior species) in two-species scenarios. As a result, you will get figures like Fig.4A in the main text.


## diversites
"beta_diversity" includes alpha, beta, and gamma diversities over the number of species in a community (2,4, ..., 10), and mean toxin sensitivity (0.1,0.2, 0.4, 0.6, 1.0). Each csv file contains alpha, beta, or gamma diversities of 100 communities (row) over the switching rate (column). From these results, we can generate Figs. 5, A3, adn A12.



## species richenss analysis
"Species_richness" contains "richness.py" and 25 csv files ("Prop_richness_speciesN_deathDelta.csv"). The python codes plots the species richness (number of persisting species) of N species communities. You can get figures like Figs.5C adn E, and Fig. A13.

"Prop_richness_speciesN_deathDelta.csv", where N is the number of species and Delta is 10x mean toxin sensitivity, contains the distributions of species richness used in the main text. Each row corresponds to each switching rate while each column represents each species richness S (S=0,1,2,...,N). 


# Appendix 2.py
This code prducces two figures (Figs. A5 and A6) in Appendix 2 (related to the scenarios in the absnece of demographic noise but the presence of the enviornmental switching).

# Apendix7.ipynb
This notebook shows how we performed the analysis in Appendix 7: (1) sub-sampling species from ten species communities to calculate the exclusion of the fittest in each pair, and (2) quantifiing the similarity between the mean probability of exclusion of the fittest and beta diversity of the whole community. 

# NoToxin
This folder contains the data where species do not interact with toxin (i.e., toxin sensitivity = 0), which is equivalent with the cases wehre we remove the toxin from our scenario. Then, one can produce Fig. A22. One will see that no toxin shows quantitatively similar results with toxin sensitivity =0.1.

# Neutral 
In this folder, we can obtain the data where two species are identical except for their label (see Fig. A21). In addition, the analysis of the neutral scenario can be used to see how the initial abundance of species affects the extinction probability in mono-culture by assuming that two species are labled as an ideintical species (Fig.A20).
