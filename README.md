# COVID-19-2D-Grid-Spread-Simulation-An-Agent-Based-Model

---

This code simulates the effect of **agent movement speeds** on the spread of COVID-19 within a **2D grid environment**. It starts with an initial population of "500 agents", randomly placed on a 10x10 grid with assigned health conditions (susceptible or infected). Agents wander the grid for a duration of "115 days", where each time step represents an hour. When an agent encounters another agent in the same grid square, there is a base probability that the susceptible agent will become infected.

Infected agents remain mobile for "6 days" before being quarantined at their current location. Post-quarantine, they are introduced back into the population as immune, with a 0.0001 probability of still spreading the virus. The simulation includes a 10x10 grid plot on the left, depicting the environment with agents as colored dots representing their health conditions, and a line graph on the right showing the proportion of infected agents over time, with a bar indicating the current population size.

To run the simulation, select the agent movement speed (v) by entering the desired speed and then click the "Run" button in the MATLAB Editor, or type `Lee_Finial_Project` in the Command Window. 

You can also modify other values, such as "N" (the initial population) or "max_infection" (the duration of infection before quarantine).
