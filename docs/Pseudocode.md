Here, we present heavily simplified pseudocode for the particular Matlab files of **MagneticFriction** to highlight the basic principles behind the simulations. However, please see the [Readme](../README.md) and the supporting information of [1] for more informations.

**CalcTorque.m** 

```pseudo

(Stokes friction torques) = - friction_coefficient x (angular velocities)

for each magnetic moment (i,j) on the slider {

    calculate the magnetic fields acting on (i,j) due to the other slider moments
    calculate the magnetic fields action on (i,j) due to the sliders on the substrate

    (magnetic torque on (i,j)) = vector_product[(magnetic moment of (i,j)), (total magnetic field acting on (i,j))]
}

(total torques) = (Stokes friction torques) + (magnetic torques)

```

**CalcVelocity.m** 

```pseudo

(angular velocities of the moments at t + time_step) = [(orientations at t + time_step) - (orientations at t)] / time_step

```

**Translation.m** 

```pseudo

(new spatial positions of the moments) = (old spatial position of the moments) + (translation velocity) x time_step

```
**Move.m** 

```pseudo

(orientations at t + time_step) = 2 x (orientation at t) - (orientation at t - time_step) + (Stokes friction torque) x time_step^2 + (magnetic torques) x time_step^2

```
