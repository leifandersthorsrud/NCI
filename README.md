# NCI
This repository contains replication codes for the research article: "Words are the new numbers: A newsy coincident index of the business cycle". 

The codes are written in MATLAB, and show how to estimate the mixed-frequency time-varying Dynamic Factor Model (DFM)
developed and used in the article referred to above.

If you use these codes, please use the following reference: 
- Thorsrud, Leif Anders, 2018. "Words are the new numbers. A newsy coincident index of the
business cycle," Journal of Business & Economic Statistics, forthcoming.

# Background: 
Decision makers and forecasters need to assess the state of the economy in real time to devise appropriate policy responses and condition on an updated information set. However, in real time, our main measure of economic activity, GDP growth, is not observed as it is compiled on a quarterly frequency and published with a considerable lag, usually up to at least one quarter.

The coincident index proposed in the article referred to above uses daily data, together with quarterly GDP growth, to provide a daily estimate of the business cycle. However, unlike other indexes with the same purpose, it has one distinctive property: The information set used to derive the index consists of daily newspaper topics.

In short, the newspaper corpus is first decomposed into distinct daily news topics. In turn, these news topics are used together with Gross Domestic Product (GDP) to derive the daily index. In the research article it is shown that the derived index has very good classification and nowcasting properties for the Norwegian business cycle. 

# Usage of the code: 
Due to the usage of proprietary data in the original research article, we can not provide data for the news topics in this repository. Instead, the repository contains codes for replicating the simulation experiment conducted in online Appendix C.1 of the research article.

The main file is called Main.m. The simulation data used for estimation can be found in the file data.mat. 
Output is saved to file. Figures similar to those reported in online Appendix C.1 can be produced by running the Graph.m script.

Please consult the license.txt file for information about redistribution, modification, and general usage of these codes.   
