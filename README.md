# abm_sir
an ABM SIR model with individual beta

The SIR model is one of the simplest compartmental models, and many models are derivatives of this basic form. 

The model consists of three compartments: 
S for the number of susceptible, I for the number of infectious, and R for the number of recovered or deceased (or immune) individuals. 
(This compartment may also be called "resistant" or "removed.")
This model is reasonably predictive[citation needed] for infectious diseases that are transmitted from human to human, 
and where recovery confers lasting resistance, such as measles, mumps and rubella.

We have seven different modality (topology)to trasmit infection:
Each agent can infect another agent (mean field)
Each agent can infect its eight (or four) immediate neighbors.
The infection is transmitted over a network, using 4 different topological structure:
random network
circle
small world
scale free network

We have different parameters for each topological structure, see list of parameters for detail
The main interface of the model is 
![Interface](./Images/abm_sir_interface.png)

##Structure 
The number of agent is N
Each agent can be in one of three states (in the course of dynamics in general, each agent can change state several times): S (susceptible), I (infected), R (recovered).
Initial conditions: no agents R, M<N agents I, the rest susceptible.
Each agent have a different capacity to infected another agent: the name of this parameter is BETA
During the setup of the model, each agent receive a different value for BETA, following a value distribution. We study the effect of different distribution of BETA. 
We call beta_distribution_schema the current BETA distribution. Actually,
we have two schema
```
if beta_distribution_schema = 1 [
  set beta_distribution [[ 0.5 1]]
  set beta_variance 0.001
 ]
 if beta_distribution_schema = 2 [
  set beta_distribution [[ 0.25 0.5][0.75 0.5]
  set beta_variance 0.01
 ]
 
 if beta_distribution_schema = 1 
    or beta_distribution_schema = 2
	[
     ask turtles [
       set beta random-normal (first rnd:weighted-one-of-list beta_distribution [ [p] -> last p ]) beta_variance
       set beta (max (list beta 0))
       set beta (min (list beta 1))
    ]
  ]
 ```
 beta_distribution_schema = 1
 we have only a value 0.5, with probabilty equal 1
 we pick BETA from a normal distribution, with 0.5 as mean and 0.001 as variance.
 
 beta_distribution_schema = 2
 we have two values, each value with probability equal 0.5
 we pick BETA from a normal distribution, with 0.5 or 0.75 as mean and 0.01 as variance.
 
##Dynamics 
divided into T_MAX time units; however, each time unit is made up of N elementary steps.

At each elementary step a random agent is extracted, let's call it k:
if k is an S,  one of its neighbors, say j, (i.e. neighbors depend from the modality: mean field, 4neighbors, network ...) is extracted,
 and if this neighbor is an I (infected), with probability p
```math
p = {\sqrt{\Beta_k * Beta_j}}
```
 k becomes I too.

 
If agent k is an I, two things are done: 
a) a neighbor is extracted, say j,(as before), if this neighbor is an S, with probability 
```math
p = {\sqrt{\Beta_k * Beta_j}}
```
also this neighbor becomes I,
b) then with probability GAMMA agent k becomes R.

Finally, if the agent k is an R, with probability RHO becomes S.




