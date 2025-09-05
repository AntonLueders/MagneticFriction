Here, we present heavily simplified pseudocode for the particular Matlab files of **MagneticFriction** to highlight the basic principles behind the simulations. However, please see the [Readme](../README.md) and the supporting information of [1] for more informations.

**CalcVelocity.m** 

```pseudo

(angular velocity of the moments at t + time_step) = [(orientations at t + time_step) - (orientations at t)] / time_step

```

**Translation.m** 

```pseudo

(new spatial positions of the moments) = (old spatial position of the moments) + (translation velocity) x time_step

```
**Move.m** 

```pseudo

(orientations at t + time_step) = 2 x (orientation at t) - (orientation at t - time_step) + (Stokes friction torque) x time_step^2 + (magnetic torques) x time_step^2

```
