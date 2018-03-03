## Motivating questions 

Plant performance often varies over space, but species often vary in how they respond to the same underlying variation. In this project I aim to describe *which* environmental layers plant performance (as meaured by the seed production of 17 annual plant species in the absence of competitors) responds to across Sedgwick reserve. This analysis will be done across the whole dataset ("Which environmental layers are annual plants responding to at Sedgwick?"). 

The second question will be based on the observation in nature (and in our dataset) that plant species don't all have the same response to environmental variation-- for example, some species achieve their maximal performance in different places. To explore this further, I will ask whether there is significant species variation in their response to the same environmental conditions (their "functional response"), and which environmental layers these species are segregating on. 

The third question will explore whether species functional traits help explain variation in species' functional response to the environment. 

----------------------------
  
## 22 Feb 2018- more thinking about the motivation

Since there's a ton of dimensionality to these models, and since I'm getting a little distracted by just running various models, I'm taking yet another step back and jotting down some thoughts on the aims, modeling approaches, potential ways to focus our analyses and simplify the data set, and limitations of the data set. To do this it's worth going back to the motivations for this study. 

## Motivation in functional trait ecology  
[Shipley et al. 2016](https://link.springer.com/article/10.1007/s00442-016-3549-x) note that a "Foundational claim" of the functional trait paradigm is that "Functional traits show general predictive relationships to measurable environmental gradients." In other words, similar environmental conditions should tend to select for a similar distribution of traits (of "response traits"). This prediction is frequently borne out by observational studies that find shifts in community-weighted mean trait values along gradients- for example, Nathan's 2008 Science paper showed turnover in functional traits within Yasuni in ways that "made sense" given the abiotic environmental variation. But on the whole, we lack generalizable relationships that underlie such patterns- for example, although 'more fertile' soils are predicted to be dominated by species with 'resource-acquisitive' traits, it's unclear how this breaks down demographically. And it's important to get this right- because ecologists have, in the past, predicted (assumed?) that community-weighted mean traits signify the adaptive value of traits at an environment, and extend this to make predictions about the performance of species along environmental gradients (e.g. "species with traits matching the CWM trait of a site will perform well at that site, if they manage to disperse in".) As I wrote in my quals proposal, this can set us up with a few competing hypotheses:  

**H0**: Species may respond identically to environmental variation, irrespective of functional traits. We may expect this pattern if the environmental gradient spans a single stress gradient axis such that all species perform poorly at some sites and well at others. In this case, observed patterns of species turnover across the landscape at Sedgwick may be driven by poor dispersal.   
**H1**: Species may have idiosyncratic responses to the environmental gradient, but the response may be unrelated to functional traits. There is potential for environmental variation to promote species coexistence, but the stability of coexistence would be uncorrelated to species' functional similarity.  
**H2**: Variation in responses of species to environmental gradients may be explained by their functional traits. This pattern might arise if the functional traits are indicative of species' environmental preference along the gradient sampled.  
**H3**: Similarity in some traits may correlate with similarity in response to certain environmental variables. For example, drought tolerance may predict demographic responses to an aridity gradient but not to a light gradient, and seed size may predict responses to a soil texture gradient but not an aridity gradient. 

![](figs/ch01-hypotheses.png)

**Note**: Having recently read parts of Hilborn and Mangel's *The Ecological Detective*, I feel compelled to speak a bit in their language. The fundamental *hypothesis* of trait based ecology is that measurable functional traits, which capture some aspect of a species' ecological strategy, influence the distrubition and abundance of the species in predictable ways. With sufficient information about a species' traits, we should be able to predict its distribution, abundance, and response to the environment; and conversely with enough information regarding the environment of a location, we should be able to predict the distribution of traits that are found in that place. There's a number of competing *models* we can use to test various parts of this framework; my goal here is to build those reasonable models and test them against one another. We have predictions about which will best describe what is going on in the environment, and model selection will either hold up our expectations or suggest alternatives.

-------------------

Before getting too deep into the specifics of the data we have from Sedgwick, I want to consider a recent paper by [Laughlin et al.](http://onlinelibrary.wiley.com/doi/10.1111/ele.12781/full) titled "Survival rates indicate that correlations between community weighted mean traits and environments can be unreliable estimates of the adaptive value of traits." Laughlin et al. address the same core questions as the ones I am interested in, by testing whether three core traits- SLA, SRL, and phenology- a) have CWM turnover across environment, and b) whether models of species survival as a function of size, environment, traits, and the trait X environment interaction yield the same ecological interpretation as the CWM results. Laughlin sets up four possible outcomes in this analysis:  

  1. There is now shift in CWM along the gradient, and there's no influence of the trait on vital rates along the gradient.  
  2. CWM shifts along the gradient, but the effect of traits on vital rates is dependent on the environment.  
  3. CWM does not shift along the gradient, but the effect of traits on vital rates appears to depend on the environmental context.  
  4. Both CWM and trait-vital rate analyses suggest adaptive value of the trait. 
  
The contrasting scenarios are summarized neatly in their Table S1:  

![](figs/laughlin-tabs1.png)

I'll come back to their methods in a bit, but for now I focus on their results and how they interpret things.

### SLA
CWM SLA was not related to sand content or to C:N ratio. In the survival analysis, the interaction term between SLA and sand content was not significant, but there was a significant interaction between SLA and C:N ratio- in High C:N ratio, survival was highest for species with low SLA and conversely. This interaction was a "strong" interaction (more on this later.)  

### SRL
CWM SRL was positively related to sand content. In survival analysis, there was a significant interaction between SRL and soil sand content, such that survival was highest for species with low SRL in soil with low sand content. This is in agreement, with the caveat that the survival analysis result was not that of a  "strong" interaction. CMW SRL was also negatively to soil C:N ratio, the interaction between SRL and soil C:N in survival analysis was not significant. 

### Flowering date
CWM flowering date was positively related to sand content. In the survival analysis, there was a significant interaction between flowering date and sand content, such that survival was higher for species with later flowering dates in sandy soils and lower for species with early flowering in sandy soil. This was a "strong" interaction. CWM flowering date was also negatively related to soil C:N ratio, but the interaction term between flowering date and soil C:N ratio was not significant, which is a conflicting result. 

### Discussion
No trait exhibited independent main effects of survival because the adaptive value of traits depended on the environmental context. 

-----------------------------

This paper is a really neat framing of part of what I'm hoping to do here- and I think I can do a neat analysis, even incorporating the ITV data I collected in 2017. Figure 1 from my dissertation proposal (see above) fits in neatly with Figure 1 of the Laughlin et al. paper- with the exception that I wasn't thinking about CWM traits during my proposal. An extention of this idea about the "adaptive value of traits" is that species "shift" their functional traits along environmental gradient in adaptive ways-- in other words, if there is predictable ITV along a gradient, then it ought to reflect the shifts in the adaptive value of the trait along the gradient. For instance, if species with high SLA respond positively to soil sand content, then we might expect individuals of a species growing in ghigher sand content to have higher SLA than individuals of the same spcies growing elsewhere. The following figure is an attempt at explaining this idea:  

![](figs/hypothesis-w-itv.png)


Of course, this figures represents the "rosiest" view of what's going on- it traits were perfectly adaptive, and if the signal came through in CWM traits, the first two panels would work out as drawn, and if ITV were perfectly adaptive, all species would display ITV as in the top panel. I suspect we are unlikely to find such a clean pattern- and more thinking will likely highlight a nuanced story. Piecing together all the possible combinations will take me some time- but one way to organize this would be to take Table S1 from Laughlin et al (as above), and in turn split the four quadrants up into three- representing the three different ITV options. 
**To Do**: Break down the quadrants in Laughlin et al. into three further possibilities of ITV. 


The nice thing about this framework is that the second two parts stand more or less independently of the first- there's interesting conclusions to be drawn from analysing the traitXenvironment drivers of vital rates, and asking whether ITV "tracks" the trait shifts you'd expect from the former. We can make a similar diagram, with "Trait X Environment effect on vital rate exists" on one axis, and "Traits shift in the same way as predicted", "Traits don't shift", and "Traits shift in opposite directions" as options on the other axis. 
**To Do**: Make a two-by-three matrix of "tXe effects exist (or not)" by "itv exists (or not), and is in the direction expected (or not)"

------------------------------   

## Building good models.

I've realized that the data we have is very non-normal, heterogeneous, may follow a negbin distribution, etc., and that my simple efforts at using sophisticated GLMMS (see `code/motivated_analyses.R` for my sandbox) are not quite sufficient. I'm going to seriously take up [Zuur's book on GLMMs for ecologists](https://link.springer.com/book/10.1007%2F978-0-387-87458-6) to get my facts right, because I think that the idea I'm developing is solid and a good analysis of it will make for a solid paper. 

The basic motivation here is:  

1. One of the fundamental assumptions of trait based functional ecology is that functional traits a) show general predictive relationships with environmental gradients because b) environments select on functional traits, such that species with the most adaptive traits are the fittest at a site.  
2. A follow-up that we still don't quite understand is whether intraspecific trait variation (ITV) is generall adaptive, or whether it is random with respect to environment (especially environment alone- usually ITV in response to environment is confounded with ITV in response to competitors.).  
3. Both of these are not very well tested. Here, we can test assumption 1a and 1b (especially 1b) by testing whether species functional traits explain how their vital rates vary across environment, , and begin to assess assumption 2 for some traits (LDMC, SLA, SRL, plant height [although how much this is a 'trait' in an ITV context is unclear to me]).


Given that I have (I think) a pretty solid conceptual basis, my main task right now is to figure out the best way to deal with the high dimensionality of the dataset, etc. So, I will (maybe "yet again") summarize here the data available to me:   

1. The demography of 17 species at 24 sites ("demography" here means i) germination rate; ii) seed production in the absence of competitors; iii) seed production in the presence of competitors)  
2. Trait data on all species (I will focus on SLA, leaf size, SRL, Max Height, Flowering Phenology, and seed mass). I have species means measured in a single site, and also have site-specific measurements of SLA, leaf size, and SRL.  
3. Environmental data from each site- the primary axes of abiotic variation are soil Ca:Mg ratio, soil nitrate concentration, and soil sand content. The biotic environment at each site is captured by an NMDS on the bacterial, fungal, and small eukaryotes. 

### Finalizing the demography dataset  

A major challenge for me right now is that the demography is a little bit all over the place. One strategy might be to focus on running the ML demography models and just work with the central tendencies of the parameter estimates- in other words, to run the ML models and end up with three tables: i) sp X site for lambda seed production; ii) sp X site for germination rate; iii) sp X site for sensitivity to competition. 

Since I ought to be pretty close to having this finalized, I will focus on that now, in the doc `ml_estimates_of_demography.Rmd`
  