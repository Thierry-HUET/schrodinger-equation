import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.linalg import solve
import time

############################### paramètres ###############################

#constante physique :
HBAR = 1
M = 1

#constante du problème :
L = 10
TEMPS = 1

#condition initiale (ici sous forme de loi normal):
SIGMAx = 1 #écart type de la fonction d'onde initial
X0 = 1.5   #espérance ou centre de la fonction d'onde
Kx0 = 1    # quantité de mouvement initial : p0 = hbar k0
SIGMAy = 1 # en y
Y0 = 1.5   # "
Ky0 = 1    # "

#potentiel periodique
TAILLE_V = 1 #taille dune barriere de potentiel
PERIODE_V = 2 #tout les combiens de distance on a une barriere
V0 = 10 #La hauteur du potentiel

#constante de numérisation :
Nx = 50 #nombre de |x>
Ny = 50
Nt = 800 #nombre de pas de temps
dx = L/Nx
dy = L/Ny
dt = TEMPS/Nt

#pour les mesures :
mesure_interval = 100

#pour l'aniamtion :
animation_interval = 5  #Intervalle pour l'animation (tout les combiens de step on sauvegarde les positions)
save_animation = False #si on sauvegarde l'animation sur la machine
save_frames = True  #si on fait une annimation

############################### fonctions ###############################

def convert(seconds):
    return time.strftime("%H:%M:%S", time.gmtime(execution_time))

def animate_wavefunction_2D(probability_for_animation):
    fig, ax = plt.subplots()
    im = ax.imshow(probability_for_animation[0], cmap="inferno", extent=[-L/2, L/2, -L/2, L/2], origin="lower")

    def update(frame_index):
        im.set_array(probability_for_animation[frame_index])
        return im,

    ani = animation.FuncAnimation(fig, update, frames=len(probability_for_animation), interval=50)
    plt.colorbar(im, label="|ψ|²")
    plt.show()

def compute_energy(psi, H): #mesure de l'energie moyenne
    return np.real(np.vdot(psi, H @ psi))  # <psi|H|psi>

def compute_U(H) : #opérateur evolution
    factor = 1j * H * dt / (2* HBAR)
    I = np.eye(Nx * Ny)
    num = I - factor
    den = I + factor
    return num, den

def potential(x, y) : # calcul le potentiel
    return (x**2 + y**2)

############################### main ###############################

start_time = time.perf_counter()

# initialisation
probability_for_animation = []

X, Y = np.meshgrid(np.linspace(-L/2, L/2, Nx), np.linspace(-L/2, L/2, Ny)) # on fait une grille 2D
# paquet d'onde
psi0 = np.exp(-0.5 * ((X - X0)**2 / SIGMAx**2 + (Y - Y0)**2 / SIGMAy**2)) * np.exp(1j * (Kx0 * X + Ky0 * Y))
psi0 /= np.linalg.norm(psi0)

#Compute Opérateur évolution
Vxy = potential(X, Y) # V(x, y)
V = np.diag(Vxy.flatten()) # matrice diagonale du potentiel dans la base des |x> tensoriel |y>

LaplacienX = (np.eye(Nx, k=1) - 2*np.eye(Nx) + np.eye(Nx, k=-1)) / dx**2
LaplacienY = (np.eye(Ny, k=1) - 2*np.eye(Ny) + np.eye(Ny, k=-1)) / dy**2
Laplacien = np.kron(LaplacienX, np.eye(Ny)) + np.kron(np.eye(Nx), LaplacienY) # np.kron = produit tensoriel

H = - (HBAR*HBAR * Laplacien / (2*M)) + V # l'Hamiltonien dans la base des |x> tensoriel |y>

U_numerateur, U_denominateur = compute_U(H) # U l'opérateur évolution = U_numerateur / U_denominateur 

# init des mesures
energy_evolution = []
mesure_time = []
"""E_moy = compute_energy(psi0, H)
print("energie initiale = ", E_moy)"""

# compute
psi = psi0.flatten()
for n in range(Nt) :
    # resoud la PDE
    psi = solve(U_numerateur, U_denominateur @ psi) # U_d|psi_n+1> = U_n|psi_n>
    
    # mesure des moyennes d'énérgie
    """E_inst = compute_energy(psi, H)
    E_moy += E_inst
    if (n % mesure_interval == 0) :
        # mesure de l'énergie
        E_moy /= mesure_interval # fait une moyenne
        energy_evolution.append(E_moy)
        E_moy = 0

        # temps ou on a fait la mesure
        mesure_time.append(n*dt)"""
        
    if n % animation_interval  == 0 : #on sauvegarde les positions pour les annimées tout les "animation_interval" pas 
        probability_for_animation.append(np.abs(psi.reshape(Nx, Ny))**2)
        print(f"image {n/animation_interval}/{Nt/animation_interval}")
    
    """normalisation = np.sum(np.abs(psi.reshape(Nx, Ny))**2)
    if np.abs(1 - normalisation) > 0.001 :
        print(f"Resultat diverge l'état n'est plus normalisé : Somme des |psi|² = {normalisation}")
        break"""

#print("energie finale = ", energy_evolution[-1])

# temps que met le programme
end_time = time.perf_counter()
execution_time = end_time - start_time
print(f"Programme exécuté en : {convert(execution_time)}")

# plot le potentiel en 3D
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Vxy, cmap='plasma')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('V(x, y)')
ax.set_title('Potentiel en 2D')
plt.show()

# lance l'animation
if save_frames and len(probability_for_animation) > 1 :
    animate_wavefunction_2D(probability_for_animation)
