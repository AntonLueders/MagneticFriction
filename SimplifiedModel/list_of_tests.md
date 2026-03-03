# Conducted Tests

In the following, we analyze the robustness of the results against changes in the numerical calculations. If one of the parameters or lines in the code is changed, everything else is kept as described in the Supplementary Information of [1] (unless stated otherwise). 

Below, **"gives the same results"** means that the average magnetic friction shows a perfect qualitative match with small quantitative deviations, and the dynamics of all three regimes are qualitatively the same (we do not claim that the degrees of freedom take exactly the same path).

## Changing the Tolerances

Decreasing the tolerances does not change the findings. The smallest tested tolerances (`options = odeset('RelTol', 1e-13, 'AbsTol', 1e-15)`) give the same results as the original analysis.

## Changing the ODE Solver

Using different ODE solvers requires different tolerances to give the same results as the original analysis. Below is a list of the tested solvers and the tolerances needed to reproduce the findings:

| Solver | Required Tolerances (`RelTol`, `AbsTol`) | Notes |
| :--- | :--- | :--- |
| `ode89` | `1e-3`, `1e-6` | Same results (tested up to `1e-13`, `1e-15`) |
| `ode113`| `1e-5`, `1e-8` | Same results |
| `ode78` | `1e-4`, `1e-7` | Same results (tested up to `1e-13`, `1e-15`) |
| `ode23` | `1e-3`, `1e-6` | Same results |
| `ode23s`| `1e-10`, `1e-13`| Same results |
| `ode23t`| `1e-7`, `1e-10` | Same results |
| `ode23tb`| `1e-9`, `1e-12` | Same results |

### Note on `ode15s`
`ode15s` does not reproduce the results of the original analysis, even for `options = odeset('RelTol', 1e-13, 'AbsTol', 1e-15)`. However, it gives qualitatively the same results: all three regimes are visible. The transition from the FM to the CP regime is shifted to $h = 8\,\text{mm}$. The magnetic friction starts at the value for the FM regime and smoothly increases throughout the entire CP regime. Similar behavior can be observed when changing the initial conditions (see below).

## Changing `tspan`

The following intervals were tested:

* `tspan = 0:0.01:10;`: Does not reproduce the findings.
* `tspan = 0:0.001:10;`: Gives the same results as the original analysis (with slightly increased magnetic friction).
* `tspan = 0:0.00005:10;`: Gives the same results.
* `tspan = 0:0.000025:10;`: Gives the same results.
* `tspan = 0:0.000015:10;`: Gives the same results.

## Changing Initial Conditions

To study the robustness of our findings against the initial conditions chosen before the equilibrium step, we scale them by a parameter $\alpha$ and check how the results change. The initial conditions are: 
$$[\varphi(0), \vartheta(0)] = \alpha \pi [1/4, 3/4]$$

Note that the equations of motion of the simplified model are not invariant under changes of the signs of the initial conditions. However, they are invariant under changing the signs of the initial conditions **and** the rotation direction of the external drive.

* **$\alpha = 0$**: Does not reproduce the findings. The system remains stuck in the FM regime (which is expected as the equations of motion are decoupled from the start).
* **$\alpha \in [10^{-7}, 1.35]$**: Reproduces the results qualitatively: All three regimes are resolved. However, the transition between the FM and CP regimes depends on $\alpha$. Decreasing $\alpha$ moves the transition towards $h = 8\,\text{mm}$. For $\alpha > 1$, the transition shifts to a larger $h$ compared to our original findings. The transition between the CP and AFM regimes is unaffected by $\alpha$.
* **$\alpha \in [-0.1, -10^{-7}]$**: Reproduces the results qualitatively. It behaves the same as for small positive $\alpha$: The transition between FM and CP occurs at $h = 8\,\text{mm}$.
* **$\alpha \in [-0.3, -0.15]$**: All regimes are present, but two distinct peaks are found: The peak at lower $h$ corresponds to the CP regime. For larger $h$, there is a peak which shows FM behavior.

## Changing the Equilibration

If the equilibration step (where the time-dependence of the external drive is excluded) is removed, but the initial conditions are changed manually to the expected results of the equilibration, the system gives the same results as the original analysis. 

We also checked if the kink of the external drive after the equilibration time affects the results. For this, we performed two distinct integrations:
1. The ODE solver is used separately for the equilibration.
2. The results are then inserted as the initial conditions for a second run corresponding to the time-dependent section of the external drive. 

This procedure gives the same results as the original analysis.

## Changing the Observation Time

We tested the stability of the dynamics in the CP regime when increasing the observation times. We tested the stability of the dynamics in the CP regime by increasing the observation times up to `tspan = 0:0.0001:100;`. Even for the largest tested time span, all three regimes still exhibit the expected dynamics and give the same results as the original analysis. However, to obtain reliable results at these large times, the tolerances must be adjusted to `options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9)`.