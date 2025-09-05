Here, we present heavily simplified pseudocode for the particular Matlab files of **MagneticFriction** to highlight the basic principles behind the simulations. However, please see the [README](../README.md) and the supporting information of [1] for more informations.



**MagneticArray.m**

```pseudo

definition and declaration of input parameters

for (slider-substrate separation) in scan_value

    initialization of the simulation (initial orientations are calculated, spatial positions of the moments are calculated, etc.)
    
    for time smaller than max_time {

        calculates torques acting on the slider moments using CalcTorque.m

        rotate the magnetic moments on the slider using Move.m

        calculate angular velocities for next simulation step using CalcVelocity.m

        if time larger than equilibration time {

            translate the slider using Translate.m
        }

        if observables should be calculated {

            if time matches (time for write output) {

                calculate observables using CalcForce.m, CalcTorqueGrid.m, CalcHysteresis.m, and CalcEnergy.m
            }    
        }
    }

}
```



**CalcTorque.m** 

```pseudo

(Stokes friction torques) = - friction_coefficient x (angular velocities)

for each magnetic moment (i,j) on the slider {

    calculates the magnetic fields acting on (i,j) due to the other slider moments

    calculates the magnetic fields action on (i,j) due to the substrate moments

    (magnetic torque on (i,j)) = vector_product[(magnetic moment of (i,j)), (total magnetic field acting on (i,j))]
}

(total torques) = (Stokes friction torques) + (magnetic torques)

```



**Move.m** 

```pseudo

(orientations at t + time_step) = 2 x (orientations at t) - (orientations at t - time_step) + (Stokes friction torques) x time_step^2 + (magnetic torques) x time_step^2

```



**CalcVelocity.m** 

```pseudo

(angular velocities of the moments at t + time_step) = [(orientations at t + time_step) - (orientations at t)] / time_step

```



**Translation.m** 

```pseudo

(new spatial positions of the moments) = (old spatial position of the moments) + (translation velocity) x time_step

```



**CalcForce.m** 

```pseudo

for each magnetic moment (i,j) on the slider {

    calculates the magnetic forces acting on (i,j) due to the other slider moments

    calculates the magnetic force action on (i,j) due the substrate moments
}

forces = sum[(total magnetic forces acting on all slider moments)]

```



**CalcTorqueGrid.m** 

```pseudo

for each magnetic moment (i,j) on the slider {

    calculates the magnetic fields action on (i,j) due to the substrate moments
    
    (magnetic torque on (i,j)) = vector_product[(magnetic moment of (i,j)), (total substrate field acting on (i,j))]
}

```



**CalcHysteresis.m** 

```pseudo

(sum of magnetic fields) = 0
(sum of torques on sublattice 1) = 0
(sum of torques on sublattice 2) = 0

for each magnetic moment (i,j) on the slider {

    calculates the magnetic fields action on (i,j) due the substrate moments

    (sum of magnetic fields) = (sum of magnetic fields) + (total substrate field acting on (i,j))

    if (i,j) is on sublattice 1 {
        (sum of torques on sublattice 1) = (sum of torques on sublattice 1) + vector_product[(magnetic moment of (i,j)), (total substrate field acting on (i,j))]
    }
    if (i,j) is on sublattice 2 {
        (sum of torques on sublattice 2) = (sum of torques on sublattice 2) + vector_product[(magnetic moment of (i,j)), (total substrate field acting on (i,j))]
    } 
}

calculates the average field per slider moment using (sum of magnetic fields)
calculates the average torque on sublattice 1 using (sum of torques on sublattice 1)
calculates the average torque on sublattice 2 using (sum of torques on sublattice 2)
calculates the magnetization of sublattice 1
calculates the magnetization of sublattice 2

```



**CalcEnergy.m** 

```pseudo

(energy of the slider) = 0
(slider-substrate interaction energy) = 0

for each magnetic moment (i,j) on the slider {

    calculates the magnetic fields acting on (i,j) due to the other slider moments
    (energy of the probe) = (energy of the probe) - dot_product[(magnetic moment of (i,j)), (magnetic field due to slider moments acting on (i,j))]

    calculates the magnetic fields action on (i,j) due the substrate moments
    (slider-substrate interaction energy) = (slider-substrate interaction energy) - dot_product[(magnetic moment of (i,j)), (magnetic field due to substrate moments acting on (i,j))]
}

```

[1] *Nonmonotonic Magnetic Friction from Collective Rotor Dynamics* by Hongri Gu, Anton Lüders, and Clemens Bechinger.
