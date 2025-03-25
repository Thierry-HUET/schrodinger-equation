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
Lx = 6
Ly = 6
TEMPS = 1
SCREEN_POSITION = -Lx/3 # proche de x = -L/2

#condition initiale (ici sous forme de loi normal):
SIGMAx = 1 #écart type de la fonction d'onde initial
X0 = 2   #espérance ou centre de la fonction d'onde
Kx0 = 15   # quantité de mouvement initial : p0 = hbar k0
SIGMAy = 1 # en y
Y0 = 0   # "
Ky0 = 0    # "

#boundary condition
DIRICHLET = False # psi au bord = 0 => potentiel infini en dehors de la boite
NEUMANN = False # d psi/dt = 0 au bord => rebond au bord
PERIODIQUE = False # la particule sort a droite => elle rentre a gauche, creer un espace infini
ABSORPTION = True #la particule sort est absorbé

# Paramètres des fentes de Young
SPLIT_CENTER1 = 0.5  # Position en y de la première fente
SPLIT_CENTER2 = -0.5  # Position en y de la deuxième fente
SLIT_WIDTH = 0.5  # Largeur des fentes
BARRIER_WIDTH = 0.3  # Épaisseur de la barrière
V0 = 20  # Hauteur de la barrière

#constante de numérisation : #Nx = Ny = 50 et Nt = 900 donne 7min de calcul #Nx = Ny = 80 et Nt = 900 donne 1h de calcul
Nx = 40 #nombre de |x>
Ny = 40 
Nt = 900 #nombre de pas de temps
dx = Lx/Nx
dy = Ly/Ny
dt = TEMPS/Nt
index_x_ecran = int((SCREEN_POSITION + (Lx/2))*(Nx/Lx))  

#pour les mesures :
mesure_interval = 100

#pour l'aniamtion :
animation_interval = 1  #Intervalle pour l'animation (tout les combiens de step on sauvegarde les positions)
save_animation = True #si on sauvegarde l'animation sur la machine
save_frames = True  #si on fait une annimation

############################### fonctions ###############################

def convert(seconds):
    return time.strftime("%H:%M:%S", time.gmtime(execution_time))

def animate_wavefunction_2D(probability_for_animation, Vxy):
    fig, ax = plt.subplots()
    extentX = [-Lx/2, Lx/2]
    extentY = [-Ly/2, Ly/2]
    extent = extentX + extentY

    # Création de l'image pour le potentiel et l'ecran en arrière-plan 
    screen_matrix = np.zeros((Ny, Nx)) 
    screen_matrix[:, index_x_ecran] = V0  # On met la valeur V0 à la colonne correspondant à l'écran pour quelque chose d'homogene
    ax.imshow(Vxy + screen_matrix, cmap="plasma", extent=extent, origin="lower", alpha=0.5)

    # Création de l'image pour l'animation
    im = ax.imshow(probability_for_animation[0], cmap="inferno", extent=extent, origin="lower", alpha=0.8)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Propagation de la fonction d'onde en 2D")

    def update(frame_index):
        im.set_array(probability_for_animation[frame_index])
        ax.set_title(f"Propagation de la fonction d'onde Frame {frame_index}/{len(probability_for_animation)-1}")
        return im,

    ani = animation.FuncAnimation(fig, update, frames=len(probability_for_animation), interval=50)
    plt.colorbar(im, label="|ψ|²")

    if save_animation:
        ani.save("young_s-interference.mp4", writer="ffmpeg", fps=20)
        plt.close(fig)
    else:
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
    V = np.zeros_like(x)

    # Définition de la barrière centrale (avec fentes)
    mask_barrier = (np.abs(x) < BARRIER_WIDTH / 2)  # Zone de la barrière en x
    mask_slit1 = np.abs(y - SPLIT_CENTER1) < SLIT_WIDTH / 2  # 1ère fente
    mask_slit2 = np.abs(y - SPLIT_CENTER2) < SLIT_WIDTH / 2  # 2ème fente

    # Ajout du potentiel
    V[mask_barrier] = V0  # Met la barrière partout
    V[mask_barrier & mask_slit1] = 0  # Fente 1
    V[mask_barrier & mask_slit2] = 0  # Fente 2

    return V

def absorbing_potential(x, y, absorption_strength=10, border_width = 0.2):
    V_abs = np.zeros_like(x)

    # Distance aux bords (normalisée entre 0 et 1 dans la zone absorbante)
    edge_x = np.minimum(np.abs(x - (-Lx / 2)), np.abs(x - (Lx / 2)))
    edge_y = np.minimum(np.abs(y - (-Ly / 2)), np.abs(y - (Ly / 2)))

    # Fonction d'absorption exponentielle sur les bords
    absorption_x = np.exp(- (edge_x / border_width) ** 2)
    absorption_y = np.exp(- (edge_y / border_width) ** 2)

    # Le potentiel absorbe dans les deux directions
    V_abs = 1j * absorption_strength  * (absorption_x + absorption_y)

    return V_abs

############################### main ###############################

start_time = time.perf_counter()

# initialisation
probability_for_animation = []
profiles_on_screen = []

X, Y = np.meshgrid(np.linspace(-Lx/2, Lx/2, Nx), np.linspace(-Ly/2, Ly/2, Ny)) # on fait une grille 2D

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

if ABSORPTION : #si les bords absorbe
    H = H.astype(complex) #car V_abs complexe, on a besoin de H complexe pour les additionnées
    V_abs = absorbing_potential(X, Y) #potentiel qui fait en sorte que la particule est absorbé aux bords
    H += np.diag(V_abs.flatten())

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
    psi_reshaped = psi.reshape(Nx, Ny)

    # Extraire la fonction d'onde sur l'écran
    psi_at_screen = psi_reshaped[:, index_x_ecran]  # On prend la colonne pour x = position de lecran
    profiles_on_screen.append(np.abs(psi_at_screen)**2) # Ajouter le profil sur l'écran à la liste

    #condition au bord :
    if DIRICHLET : # Boite quantique, pour un puits de potentiel infini
        psi_reshaped[0, :] = psi_reshaped[-1, :] = 0  # en x
        psi_reshaped[:, 0] = psi_reshaped[:, -1] = 0  # en y
    elif NEUMANN : #La particule rebondit au bord sans perdre d’énergie.
        psi_reshaped[0, :] = psi_reshaped[1, :]  # Bord gauche x
        psi_reshaped[-1, :] = psi_reshaped[-2, :]  # Bord droit x
        psi_reshaped[:, 0] = psi_reshaped[:, 1]  # Bord bas y
        psi_reshaped[:, -1] = psi_reshaped[:, -2]  # Bord haut y
    elif PERIODIQUE : # Une particule qui sort d’un bord réapparaît de l’autre côté.
        psi_reshaped[0, :] = psi_reshaped[-2, :]  # Bord gauche = Bord droit x
        psi_reshaped[-1, :] = psi_reshaped[1, :]  # Bord droit = Bord gauche x
        psi_reshaped[:, 0] = psi_reshaped[:, -2]  # Bord bas = Bord haut y
        psi_reshaped[:, -1] = psi_reshaped[:, 1]  # Bord haut = Bord bas y
    psi = psi_reshaped.flatten()

    # mesure des moyennes d'énérgie
    """E_inst = compute_energy(psi.reshape(Nx, Ny), H)
    E_moy += E_inst
    if (n % mesure_interval == 0) :
        # mesure de l'énergie
        E_moy /= mesure_interval # fait une moyenne
        energy_evolution.append(E_moy)
        E_moy = 0

        # temps ou on a fait la mesure
        mesure_time.append(n*dt)"""
        
    if n % animation_interval  == 0 : #on sauvegarde les positions pour les annimées tout les "animation_interval" pas 
        probability_for_animation.append(np.abs(psi_reshaped)**2)
        print(f"image {n/animation_interval}/{Nt/animation_interval}")

#print("energie finale = ", energy_evolution[-1])
print("compute finis")

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
    animate_wavefunction_2D(probability_for_animation, Vxy)

# temps que met le programme
end_time = time.perf_counter()
execution_time = end_time - start_time
print(f"Programme exécuté en : {convert(execution_time)}")

# Affichage des bandes d'interférence
final_profile = np.array(profiles_on_screen)  # Tableau de forme (Nt, Ny)

fig, ax = plt.subplots(figsize=(8, 6))
ax.imshow(final_profile, cmap="inferno", aspect="auto", extent=[-Ly/2, Ly/2, 0, Nt * dt])
ax.set_xlabel("y")
ax.set_ylabel("Temps")
ax.set_title(f"Evolution des bandes d'interférence sur l'écran en x = {SCREEN_POSITION}")
plt.show()

# Tracer la somme des probabilités sur l'écran
plt.plot(np.linspace(-Ly/2,Ly/2, Nx), np.sum(final_profile, axis=0), color = 'r')
plt.title(f"Bandes d'interférence sur l'écran en x = {SCREEN_POSITION}")
plt.xlabel("Position en y") 
plt.ylabel("Somme des probabilités au cours du temps")
plt.show()