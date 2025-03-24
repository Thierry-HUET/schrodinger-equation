# Schrodinger Equation - Quantum Simulations  

This repository contains numerical simulations of quantum mechanics problems using **Python**. The simulations are based on solving the **time-dependent Schrödinger equation (TDSE)** for different potential configurations.  

The project is structured as an **introduction to Partial Differential Equations (PDEs)** applied to quantum physics, including examples such as the **heat equation, harmonic oscillator, quantum tunneling, and Young's double-slit experiment**.  

---

## 📌 **Project Overview**  

### 1️⃣ **Introduction: Heat Equation (PDE basics)**  
- Implements the **heat diffusion equation**, which serves as an introduction to PDE solvers.  
- Uses **finite difference methods** to evolve the system in time.  
- **Visualized with Matplotlib animations**.  

### 2️⃣ **Quantum Harmonic Oscillator**  
- Solves the Schrödinger equation for a **harmonic potential** \( V(x) = \frac{1}{2} m \omega^2 x^2 \).  
- Shows **wavefunction evolution**.  
- Implements numerical evolution using the **Crank-Nicolson method** for stability.  
- **Visualized with Matplotlib animations**.  

### 3️⃣ **Quantum Tunneling Effect**  
- Simulates a **wave packet propagating through a periodic potential**.  
- Demonstrates **barrier penetration**, a key quantum phenomenon.  
- **Visualized with Matplotlib animations**.  

### 4️⃣ **Young’s Double-Slit Experiment**  
- Simulates **wave interference in a 2D potential**.  
- Models the **probability density** behind two slits.  
- Helps visualize quantum superposition and the **wave-particle duality**.  
- **Visualized with Matplotlib animations**.  

---

## 🛠 **Installation & Dependencies**  
This project requires Python **3.7+** and the following libraries:  
```bash
pip install numpy matplotlib scipy