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
TEMPS = 100

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
Nx = 200 #nombre de |x>
Ny = 200
Nt = 20000 #nombre de pas de temps
dx = L/Nx
dy = L/Ny
dt = TEMPS/Nt

#pour les mesures :
mesure_interval = 100

#pour l'aniamtion :
animation_interval = 20  #Intervalle pour l'animation (tout les combiens de step on sauvegarde les positions)
save_animation = False #si on sauvegarde l'animation sur la machine
save_frames = True  #si on fait une annimation

############################### fonctions ###############################

def convert(seconds):
    return time.strftime("%H:%M:%S", time.gmtime(execution_time))

def animate_trajectory(probability_for_animation,energy_evolution, Vx, L): #Animation des trajectoires de particule
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    scat = ax2.scatter([], [], s=20, color ='blue')
    ax1.plot(x, Vx, color = 'g', label="V(x)",  linestyle = '--')
    ax1.plot(x, energy_evolution, color = 'red', label="<H>")
    
    ax1.set_xlim(-L/2, L/2)
    ax2.set_ylim(0, 1.1)
    
    ax1.set_xlabel('x')
    ax1.set_ylabel("V(x) & <H>", color = 'green')
    ax1.set_ylabel("<H>", color = 'red')
    ax2.set_ylabel('|ψ|²', color='blue')
    
    title = "Tunneling effect"
    ax1.set_title(title)

    def init():
        scat.set_offsets(np.empty((0, 2)))
        return (scat,)
    
    def update(frame_index):
        coords = probability_for_animation[frame_index]
        scat.set_offsets(coords)
        ax1.set_title(f"{title} Frame {frame_index}/{len(probability_for_animation)-1}")
        return (scat,)
    
    #print("Shape for animation:", probability_for_animation.shape)
    
    ani = animation.FuncAnimation(
        fig, update,
        frames=len(probability_for_animation),
        init_func=init,
        blit=True,
        interval=50 #50ms -> 20fps
    )
    if save_animation :
        ani.save("tunneling-effect-2D.mp4", writer="ffmpeg", fps=20) #pour sauvegarder l'animation en .mp4 avec ffmpeg
        plt.close(fig)
    else :
        plt.show()

def compute_energy(psi, H): #mesure de l'energie moyenne
    return np.real(np.vdot(psi, H @ psi))  # <psi|H|psi>

def compute_U(H) : #opérateur evolution
    factor = 1j * H * dt / (2* HBAR)
    I = np.eye(Nx)
    num = I - factor
    den = I + factor
    return num, den

def potential(x, y) : # calcul le potentiel
    return (x**2 + y**2)

############################### main ###############################

start_time = time.perf_counter()

# initialisation
probability_for_animation = []

x, y = np.meshgrid(np.linspace(-L/2, L/2, Nx), np.linspace(-L/2, L/2, Ny)) # on fait une grille 2D

if True : # paquet d'onde
    psi0 = [np.exp(-0.5 * (x - X0)**2 / SIGMAx**2) * np.exp(1j * Kx0 * x), np.exp(-0.5 * (y - Y0)**2 / SIGMAy**2) * np.exp(1j * Ky0 * y)]
    psi0 /= np.linalg.norm(psi0)
else : # un etat propre donc un |x> tensoriel |y> avec x = X0 et y = YO
    psi0 = np.zeros(Nx,Ny)
    psi0[int((X0+(L/2))*(Nx/L)), int((Y0+(L/2))*(Ny/L))] = 1
#Compute Opérateur évolution
Vxy = potential(x, y) # V(x, y)
V = np.diag(Vxy.flatten()) # matrice diagonale du potentiel dans la base des |x> tensoriel |y>

LaplacienX = (np.eye(Nx, k=1) - 2*np.eye(Nx) + np.eye(Nx, k=-1)) / dx**2
LaplacienY = (np.eye(Ny, k=1) - 2*np.eye(Ny) + np.eye(Ny, k=-1)) / dy**2
Laplacien = np.kron(LaplacienX, np.eye(Ny)) + np.kron(np.eye(Nx), LaplacienY) # np.kron = produit tensoriel

H = - HBAR*HBAR * Laplacien / (2*M) + V # l'Hamiltonien dans la base des |x> tensoriel |y>

U_numerateur, U_denominateur = compute_U(H) # U l'opérateur évolution = U_numerateur / U_denominateur 

# init des mesures
energy_evolution = []
mesure_time = []
E_moy = compute_energy(psi0, H)
print("energie initiale = ", E_moy)

# compute
psi = psi0
for n in range(Nt) :
    # resoud la PDE
    psi = solve(U_numerateur, U_denominateur @ psi) # U_d|psi_n+1> = U_n|psi_n>
    
    # mesure des moyennes d'énérgie
    E_inst = compute_energy(psi, H)
    E_moy += E_inst
    if (n % mesure_interval == 0) :
        # mesure de l'énergie
        E_moy /= mesure_interval # fait une moyenne
        energy_evolution.append(E_moy)
        E_moy = 0

        # temps ou on a fait la mesure
        mesure_time.append(n*dt)
        
    if n % animation_interval  == 0 : #on sauvegarde les positions pour les annimées tout les "animation_interval" pas 
        positions = []
        for i in range(Nx) :
            positions.append([i*dx - L/2, np.abs(psi[i])**2])
        probability_for_animation.append(positions)
        #print(f"image {n/animation_interval}/{Nt/animation_interval}")
    
    normalisation = np.sum(np.abs(psi)**2)
    if np.abs(1 - normalisation) > 0.001 :
        print(f"Resultat diverge l'état n'est plus normalisé : Somme des |psi|² = {normalisation}")
        break

print("energie finale = ", energy_evolution[-1])

# temps que met le programme
end_time = time.perf_counter()
execution_time = end_time - start_time
print(f"Programme exécuté en : {convert(execution_time)}")

# lance l'animation
if save_frames and len(probability_for_animation) > 1 :
    animate_trajectory(probability_for_animation,[np.mean(energy_evolution[1:])]*Nx, Vx, L)
