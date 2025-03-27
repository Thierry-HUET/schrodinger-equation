import numpy as np
import matplotlib.pyplot as plt

def diffraction_with_unique_particle(Ny, Nx, Ly, Lx, dt, wavelength, slit_distance, screen_distance, t_start=0):
    # Positions sur l'écran
    y = np.linspace(-Ly/2, Ly/2, Ny)
    
    # Intensité théorique (fonction de diffraction)
    I = -(np.cos(np.pi * slit_distance * y / (wavelength * screen_distance)))**2
    
    # Facteur d'atténuation exponentielle pour la décroissance de l'intensité
    intensity_decay = np.exp(-y**2 / 2)  # Une décroissance gaussienne
    I = I * intensity_decay  # Appliquer la décroissance
    
    # Profondeur de la probabilité (la particule commence à t_start, puis se propage)
    t = np.linspace(0, Nx * dt, Nx)  # Temps fictif
    profile = np.outer(-np.abs(t-t_start-1), I)  # Décroissance progressive des bandes
    
    return profile

# Paramètres de la simulation
Ny, Nx = 200, 100  # Résolution en y et en temps
Ly, Lx = 6, 6  # Taille de l'écran
dt = 0.01  # Pas de temps
wavelength = 0.15  # Longueur d'onde simulée
slit_distance = 2  # Distance entre les fentes
screen_distance = 10 # Distance de la fente à l'écran

# Génération du profil théorique avec particule unique
profile_with_particle = diffraction_with_unique_particle(Ny, Nx, Ly, Lx, dt, wavelength, slit_distance, screen_distance)

# Affichage
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# 1. Affichage de l'évolution temporelle des bandes d'interférence
im1 = axes[0].imshow(profile_with_particle, cmap="inferno", aspect="auto", extent=[-Ly/2, Ly/2, 0, Nx * dt])
axes[0].set_xlabel("y")
axes[0].set_ylabel("Temps")
axes[0].set_title("Évolution des bandes d'interférence (Particule Unique)")
fig.colorbar(im1, ax=axes[0], label="|ψ|²")

# 2. Tracé de la somme des probabilités sur l'écran
axes[1].plot(np.linspace(-Ly/2, Ly/2, Ny), np.sum(profile_with_particle, axis=0), color='r')
axes[1].set_title("Somme des probabilités (Particule Unique)")
axes[1].set_xlabel("Position en y")
axes[1].set_ylabel("Intensité cumulée")
axes[1].grid(True)

# 3. Projection 2D de la somme des probabilités
screen = np.tile(np.sum(profile_with_particle, axis=0), (Ny, 1))  # Réplication correcte
im3 = axes[2].imshow(screen, cmap="inferno", aspect="auto", extent=[-Ly/2, Ly/2, 0, Ny])
axes[2].set_xlabel("y")
axes[2].set_ylabel("z")
axes[2].set_title("Projection des bandes d'interférence")
fig.colorbar(im3, ax=axes[2], label="Intensité")

plt.tight_layout()
plt.show()
