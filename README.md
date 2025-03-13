## Overview

This project simulates the spread of infection using an integro-differential model within a Susceptible-Infected-Recovered (SIR) framework. The model incorporates a spatial interaction kernel to account for non-homogeneous mixing between susceptible and infected individuals in one and two-dimensional spaces. 
The project explores how infections spread across a population within a confined space, considering factors like infection rate, recovery rate, and spatial interactions. 
The model is space and time-dependent, allowing for the simulation of various conditions, including the implementation of control functions that mimic interventions like lockdowns.

## Key Concepts

*   **SIR Model:** A fundamental epidemiological model categorizing individuals into Susceptible (S), Infected (I), and Recovered (R) compartments to study disease transmission dynamics.
*   **Integro-Differential Equations (IDE):** Equations involving both derivatives and integrals, used here to model systems where dynamics are influenced by historical states and spatial interactions.
*   **Spatial Interaction Kernel:** A function, denoted as k(t, x-y), that describes the interaction between spatial points, influencing the spread of infection based on distance and control measures.
*   **Discretization:** The process of transforming continuous equations into discrete forms suitable for numerical computation, using methods like the Euler method.

## Model Description

The model is based on the SIR model, dividing the population into three categories:

*   **Susceptible (S):** Individuals prone to infection upon contact with infected individuals.
*   **Infected (I):** Individuals currently infected and capable of transmitting the disease.
*   **Recovered (R):** Individuals who have recovered and are no longer infectious.

The model uses the following key equations:

*   **Infection dynamic equation:**

    *   In the basic SIR model, the change in infected individuals $\( z \)$ over time is defined as $\(\frac{d}{dt} z(t,x) = \beta (1-z-r) z - \gamma z\)$.
*   **Recovery dynamic equation:**

    *   The change in recovered individuals $\( r \)$ over time is defined as $\(\frac{d}{dt} r(t,x) = \gamma z\)$
*   **Integro-Differential Equation for Infection (IDE):**

    *   $\(\frac{d}{dt} z(t,x) = \beta (1-z-r) \int_{0}^{1} z(t,y)k(t,x,y) \,dy - \gamma z\)$
*   **Kernel Function:**
        $K = (1 -u(t))(c*e^{-\delta|x - y|})  $
## Numerical Approximation

The integro-differential equations are solved numerically using the Euler method [13-16]. This involves discretizing the spatial and temporal domains [17, 18].

*   **Discretization:** Approximating continuous derivatives with discrete differences using a time step τ.
*   **Euler Method:** An iterative approach to approximate the solution of the differential equations at discrete time steps.
*   **One-Dimensional Discretization:**
    *   $\(z_{i}^{k+1} = z_{i}^{k} + \tau [ \beta (1-z_{i}^{k}-r_{i}^{k}) \int_{0}^{1} z(t,y)k(t,x,y) dy - \gamma z_{i}^{k} ]\)$
    *   $\(r_{i}^{k+1} = r_{i}^{k} + \tau [\gamma z_{i}^{k}]\)$
*   **Two-Dimensional Discretization:**
    *   $\(z_{i,j}^{k+1} = z_{i,j}^{k} + \tau [ \beta (1-z_{i,j}^{k}-r_{i,j}^{k}) \int_{0}^{1} \int_{0}^{1} z(t,x,y)k(t,x,y) dy dx - \gamma z_{i,j}^{k} ]\)$
    *    $\(r_{i,j}^{k+1} = r_{i,j}^{k} + \tau [\gamma z_{i,j}^{k}]\)$

## Implementation

The models are implemented using Python with the following libraries:

*   **numpy:** For numerical computations.
*   **matplotlib:** For creating visualizations and animations.


### Parameters

Key parameters that can be adjusted in the models include:

*   **β (beta):** Infection rate.
*   **γ (gamma):** Recovery rate.
*   **τ (tau):** Time step.
*   **Kernel Parameters:** Parameters within the spatial interaction kernel, such as $\( c \)$ and $\( \delta \)$, which influence the range and strength of spatial interactions.

## Results

The project provides visualizations of infection spread in both one and two-dimensional domains.  The figures generated shows the distribution of infections (z) and recoveries (r) over time and space.

*   **One-Dimensional Spatial Domain:** Assuming infection starts between points 50 and 80 in a 200-point domain.
*   **Two-Dimensional Spatial Domain:**
    *   **Case 1:** Initial infection at the center of a 30x30 grid.
    *   **Case 2:** Initial infection at the corners of the domain.

## Usage

1.  **Clone the repository.**
2.  **Install the required libraries:**

    ```bash
    pip install numpy matplotlib
    ```
3.  **Run the Python script** to simulate the infection spread and generate visualizations.
    *   `python "IDE-Python Code".py`
   
4.  **Modify the parameters** within the scripts to explore different scenarios and conditions.
