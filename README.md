# abm_sir
an ABM SIR model with individual beta

The SIR model is one of the simplest compartmental models, and many models are derivatives of this basic form. 
The model consists of three compartments: 
S for the number of susceptible, I for the number of infectious, and R for the number of recovered or deceased (or immune) individuals. 
(This compartment may also be called "resistant" or "removed.")
This model is reasonably predictive[citation needed] for infectious diseases that are transmitted from human to human, 
and where recovery confers lasting resistance, such as measles, mumps and rubella.

We have two different modality
Each cell can infect its eight (or four) immediate neighbors.

One type of agent.
Each agent can be in one of three states (in the course of dynamics in general, each agent can change state several times): S (susceptible), I (infected), R (recovered).

Initial conditions: no agents R, M <N agents I, the rest susceptible.

Dynamics: divided into T_MAX time units (days, months, whatever you want); however, each time unit is made up of N elementary steps.

At each elementary step a random agent is extracted, let's call it k: if k is an S, 
one of its neighbors is extracted (i.e. any other agent if we are in the mean field, one of its 4 or 8 first neighbors in the other two),
 and if this neighbor is an I (infected), with probability BETA (see point 7) k becomes I too.

 update_infection - -

For the moment BETA is a fixed parameter of the model, but when we refine the model it will vary, over time and from agent to agent (how we will see it). If agent k is an I, two things are done: a) a neighbor is extracted (similarly to point 5), if this neighbor is an S, with probability BETA also this neighbor becomes I, b) then with probability GAMMA agent k becomes R.

Finally, if the agent k is an R, with probability RHO becomes S.

For clarity, I repeat that points 6-10 represent a single elementary step, N elementary steps make up a time unit, the whole simulation lasts T_MAX time units.

In summary: the parameters that the user can choose are BETA, GAMMA, RHO (which being probabilities go from 0 to 1, I would do it in steps of 0.05 or smaller); topology (mean field, 2D with 4 p.v., 2D with 8 p.v.) the initial M number of infected is also chosen by the user. N, on the other hand, is fixed, I will make 100 or 200 (you see, in general the more I am, the better, clearly if too many, then the dynamics are slow).
