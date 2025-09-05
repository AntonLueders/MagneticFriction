Here, we present simplified pseudocode for the particular Matlab files of **MagneticFriction** to highlight the basic principles behind the simulations. However, please see the [Readme](../README.md) and the supporting information of [1] for more informations.

**CalcVelocity.m** 

```pseudo

(Angular velocity of the moments at t + time_step) = [(Angular velocity of the moments at t + time_step) - [(Angular velocity of the moments at t)] / time_step

```

**Translation.m** 

```pseudo

(New spatial positions of the moments) = (old spatial position of the moments) + (translation velocity) x time_step

```
