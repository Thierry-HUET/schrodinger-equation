import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.linalg import solve

############################### paramètres ###############################

#constante physique :
HBAR = 1
M = 1

#constante du problème :
L = 10
TEMPS = 10

#condition initiale (ici sous forme de loi normal):
SIGMA = 1
X0 = 2

#constante de numérisation :
Nx = 1000 #nombre de |x>
Nt = 20000 #nombre de pas de temps
dx = L/Nx
dt = TEMPS/Nt

#pour l'aniamtion :
animation_interval = 25  #Intervalle pour l'animation (tout les combiens de step on sauvegarde les positions)
save_animation = False #si on sauvegarde l'animation sur la machine
save_frames = True  #si on fait une annimation

############################### fonctions ###############################

def animate_trajectory(positions_for_animation, L): #Animation des trajectoires de particule
    fig, ax = plt.subplots(figsize=(6,6))
    scat = ax.scatter([], [], s=20) #, color='blue'
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    title = "Harmonic Oscillator"
    ax.set_title(title)
    plt.ylim([0,1])
    def init():
        scat.set_offsets(np.empty((0, 2)))
        return (scat,)
    
    def update(frame_index):
        coords = positions_for_animation[frame_index]
        scat.set_offsets(coords)
        ax.set_title(f"Frame {frame_index}/{len(positions_for_animation)-1}")
        return (scat,)
    
    #print("Shape for animation:", positions_for_animation.shape)
    
    ani = animation.FuncAnimation(
        fig, update,
        frames=len(positions_for_animation),
        init_func=init,
        blit=True,
        interval=50 #50ms -> 20fps
    )
    if save_animation :
        ani.save("harmonic-oscillator.mp4", writer="ffmpeg", fps=20) #pour sauvegarder l'animation en .mp4 avec ffmpeg
        plt.close(fig)
    else :
        plt.show()

def compute_U(H) :
    factor = 1j * H * dt / (2* HBAR)
    I = np.eye(Nx)
    num = I - factor
    den = I + factor
    return num, den

#potentiel
def potential(x) :
    return x * x * x

############################### main ###############################

#initialisation
positions_for_animation = []

x = np.arange(-L/2, L/2, dx)
psi0 = np.exp(-0.5 * (x - X0)**2 / SIGMA**2)*(1/(SIGMA*np.sqrt(2*np.pi))) #np.zeros(Nx)
#psi0[Nx//4] = 1

#Compute Opérateur évolution
Vx = potential(x) #V(x)
V = np.diag(Vx) #matrice diagonale du potentiel dans la base des |x>

Laplacien = (np.eye(Nx,k=1) - 2*np.eye(Nx) + np.eye(Nx,k = -1))/(dx*dx) 

H = - HBAR*HBAR * Laplacien / (2*M) + V # l'Hamiltonien dans la base des |x>

U_numerateur, U_denominateur = compute_U(H) # U l'opérateur évolution = U_numerateur / U_denominateur 

#compute
psi = psi0
for n in range(Nt-1) :
    #resoud la PDE
    psi = solve(U_numerateur,U_denominateur @ psi) #np.matmul(U, psi)
    
    if n % animation_interval  == 0 : #on sauvegarde les positions pour les annimées tout les "animation_interval" pas 
        positions = []
        for i in range(Nx) :
            positions.append([i*dx, np.abs(psi[i])**2])
        positions_for_animation.append(positions)

        normalisation = np.sum(np.abs(psi)**2)
        if 1 - normalisation > 0.01 :
            print(f"Resultat diverge l'état n'est plus normalisé : Somme des |psi|² = {normalisation}")
            break

# lance l'animation
if save_frames and len(positions_for_animation) > 1 :
    animate_trajectory(positions_for_animation, L)