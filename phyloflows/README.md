**phyloflows** provides tools to infer transmission flows within and
between population groups from phylogenetically observed flow counts.

**phyloflows** was primarily written to estimate transmission flows from
counts of phylogenetically reconstructed likely transmission pairs, in
whom the direction of transmission could be inferred with high support.
Such input data can be obtained with
[phyloscanner](https://github.com/BDI-pathogens/phyloscanner).

Different types of input data can also be used, so long as they
represent observed counts of flow events from and to population
sub-groups. Methodologically, in **phyloflows**, there is no limitation
to deep-sequence phylogenetic data.

The primary feature of **phyloflows** is that the underlying statistical
model allows adjusting for sampling heterogeneity across population
groups while estimating transmission flows, a common problem in
molecular epidemiologic analyses.

Installation
------------

Install the latest development version from **GitHub**:

    install_github("BDI-pathogens/phyloscanner/phyloflows", dependencies=TRUE, build_vignettes=FALSE)
    require(phyloflows)

If you have issues with installation/running of phyloflows, please raise
an Issue on the GitHub page and we will get back to you.

To see beautiful latex equations in the markdown below, you may want to
[install the MathJax plugin for
Github](https://github.com/orsharir/github-mathjax).

Getting started
---------------

If you are just getting started with **phyloflows**, please take a look
at the sections below, and the tutorial vignettes at the bottom of this
page.

Our job
-------

Suppose the true number of transmission flows between two population
groups “1” and “2” is
(*z*<sub>11</sub> = 50, *z*<sub>12</sub> = 10, *z*<sub>21</sub> = 20, *z*<sub>22</sub> = 20),
 so that the true proportion of transmission flows within and between
groups is
(*π*<sub>11</sub> = 0.5, *π*<sub>12</sub> = 0.1, *π*<sub>21</sub> = 0.2, *π*<sub>22</sub> = 0.2).
 **The primary aim** is to estimate the vector of proportions
*π* = (*π*<sub>11</sub>, *π*<sub>12</sub>, *π*<sub>21</sub>, *π*<sub>22</sub>).
This looks simple. **However in real life one main complication is**
that only a proportion of the true transmission flows are observed (or
sampled), and that the sampling probabilities differ. In our experience
this situation not only occurs frequently in observational epidemiologic
studies but also controlled situations (trials). Suppose that population
“1” is sampled with probability *s*<sub>1</sub> = 0.6 and population “2”
with probability *s*<sub>2</sub> = 1. On expectation, we then observe
the transmission flows
*n*<sub>*i**j*</sub> = *z*<sub>*i**j*</sub> \* *s*<sub>*i*</sub> \* *s*<sub>*j*</sub>,
if we assume that sampling within and between the two groups is
independent of each other, and if we ignore the fact that in reality
sampling is without replacement. It won t matter much in reasonably
large populations anyway, considering the extent of typical sampling
differences. So what we observe is on expectation
(*n*<sub>11</sub> = 18, *n*<sub>12</sub> = 6, *n*<sub>21</sub> = 12, *n*<sub>22</sub> = 20),
 If we don’t adjust for the sampling, a crude estimate of the true
proportion of transmission flows would be
(*π̂*<sub>11</sub> = 18/56 = 0.321, *π̂*<sub>12</sub> = 0.107, *π̂*<sub>21</sub> = 0.214, *π̂*<sub>22</sub> = 0.357)
 and the worst case error is 50%-32.1%=17.9%. That is pretty bad. **The
job of phyloflow is to return better estimates**, with a worst case
error lower than 2-3%.

Our solution
------------

One very quick solution to the sampling problem would be to infer the
actual transmission flows from what we know, i.e.
*ẑ*<sub>*i**j*</sub> = *n*<sub>*i**j*</sub>/(*s*<sub>*i*</sub> \* *s*<sub>*j*</sub>),
 and then to calculate
*π̂*<sub>*i**j*</sub> = *ẑ*<sub>*i**j*</sub>/∑<sub>*k**l*</sub>*ẑ*<sub>*k**l*</sub>.
 But this does not address the further problems that the sampling groups
may be different to the population groups between whom we want to infer
transmission flows; that the sampling probabilities are usually not
known themselves; and that we also want interpretable uncertainty
estimates for *π̂*. So we use a Bayesian approach, and write down the
likelihood of the complete data,
*p*(*z*|*Z*, *π*) = *M**u**l**t**i**n**o**m**i**a**l*(*z*; *Z*, *π*),
 where *Z* is the total number of transmissions,
*Z* = ∑<sub>*k**l*</sub>*z*<sub>*k**l*</sub>. Next, we add the
likelihood of the sampling process,
*p*(*n*<sub>*i**j*</sub>|*z*<sub>*i**j*</sub>, *s*<sub>*i*</sub>, *s*<sub>*j*</sub>) = *B**i**n**o**m**i**a**l*(*n*<sub>*i**j*</sub>; *z*<sub>*i**j*</sub>, *s*<sub>*i*</sub> \* *s*<sub>*j*</sub>).
 Finally, we complete the model with prior distributions on the
proportions, *p*(*π*), which we generally choose in an uninformative
way; prior distributions on the sampling probabilities,
*p*(*s*<sub>*i*</sub>), which we generally choose based on other
availabe information; and a prior distribution on the total number of
transmissions, *p*(*Z*), which we also choose based on available
information. And then we use Bayes Theorem.

This looks pretty straightforward. The issue is that in real-world data
analyses, the total number of unknown parameters is 1,000 or even more.
And so we also need a nifty computational inference engine.
**phyloflow** uses a special type of Markov Chain Markov Carlo
algorithms for this job. The main feature is that it has zero tuning
variables, so the blame is 100% on us if it does not work. Just get in
touch then.

Examples
--------

1.  [Simulating data.](vignettes/01_simulating_data.md)

2.  [A very first example. Always start with
    me.](vignettes/02_basic_example.md)

3.  [Performing MCMC diagnostic checks.](vignettes/03_diagnostics.md)

4.  [Aggregating MCMC output to a new set of
    variables.](vignettes/04_aggregating.md)

5.  [Calculating sources, onward transmissions and flow
    ratios.](vignettes/05_keyquantities.md)

6.  [Testing **phyloflows** multi-level model
    and MCMC.](vignettes/06_test_sampling_adjustments.md)

7.  [Age analysis example.](vignettes/07_age_analysis.md)

8.  [Real-world example.](vignettes/08_practical_example.md)
