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
Lx = 10
Ly = 10
TEMPS = 1
SCREEN_POSITION = np.array([Lx/3 - Lx/2, Lx/4 - Lx/2, Lx/6 - Lx/2]) # proche de x = -L/2

#condition initiale (ici sous forme de loi normal):
SIGMA_X = 1 # écart type de la fonction d'onde initial
X0 = Lx/4    # espérance ou centre de la fonction d'onde
Kx0 = 15    # quantité de mouvement initial : p0 = hbar k0
SIGMA_Y = 1 # en y
Y0 = 0      # "
Ky0 = 0     # "

#boundary condition
DIRICHLET = False # psi au bord = 0 => potentiel infini en dehors de la boite
NEUMANN = False # d psi/dt = 0 au bord => rebond au bord
PERIODIQUE = False # la particule sort a droite => elle rentre a gauche, creer un espace infini
ABSORPTION = True #la particule sort est absorbé

# Paramètres des fentes de Young
X_SLIT_POSITION = Lx/8 # ou se trouve les fentes en x
SPLIT_CENTER1 = 0.4  # Position en y de la première fente
SPLIT_CENTER2 = - SPLIT_CENTER1  # Position en y de la deuxième fente
SLIT_WIDTH = 0.4  # Largeur des fentes
BARRIER_WIDTH = 0.3  # Épaisseur de la barrière
V0 = 150  # Hauteur de la barrière

#constante de numérisation : #Nx = Ny = 50 et Nt = 900 donne 7min de calcul #Nx = Ny = 80 et Nt = 900 donne 1h de calcul
Nx = 125 #nombre de |x>
Ny = 125 
Nt = 1000 #nombre de pas de temps
dx = Lx/Nx
dy = Ly/Ny
dt = TEMPS/Nt
screen_index_x = ((SCREEN_POSITION + (Lx/2)) * (Nx/Lx)).astype(int) # Indice de l'écran dans la grille 

#pour les mesures :
mesure_interval = 10

#pour l'aniamtion :
animation_interval = 5  #Intervalle pour l'animation (tout les combiens de step on sauvegarde les positions)
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
    for i in screen_index_x :
        screen_matrix[:, i] = V0/2  # On met la valeur V0/2 à la colonne correspondant aux écrans pour quelque chose d'homogene
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
    mask_barrier = (np.abs(x - X_SLIT_POSITION) < BARRIER_WIDTH / 2)  # Zone de la barrière en x
    mask_slit1 = np.abs(y - SPLIT_CENTER1) < SLIT_WIDTH / 2  # 1ère fente
    mask_slit2 = np.abs(y - SPLIT_CENTER2) < SLIT_WIDTH / 2  # 2ème fente

    # Ajout du potentiel
    V[mask_barrier] = V0  # Met la barrière partout
    V[mask_barrier & mask_slit1] = 0  # Fente 1
    V[mask_barrier & mask_slit2] = 0  # Fente 2

    return V

def absorbing_potential(x, y, absorption_strength=15, border_width = 0.75):
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
profiles_on_screen = [[] for _ in range(len(SCREEN_POSITION))]

X, Y = np.meshgrid(np.linspace(-Lx/2, Lx/2, Nx), np.linspace(-Ly/2, Ly/2, Ny)) # on fait une grille 2D

# Initialisation du paquet d'onde
psi0 = np.exp(-0.5 * ((X - X0)**2 / SIGMA_X**2 + (Y - Y0)**2 / SIGMA_Y**2)) * np.exp(1j * (Kx0 * X + Ky0 * Y))
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
normalisation_evolution = []
mesure_time = []
E_moy = compute_energy(psi0.flatten(), H)
normalisation_moy = np.sum(np.abs(psi0.flatten())**2)
print("energie initiale = ", E_moy, "la norme initiale = ", normalisation_moy)

# Propagation temporelle
psi = psi0.flatten()
for n in range(Nt) :
    # resoud la PDE
    psi = solve(U_numerateur, U_denominateur @ psi) # U_d|psi_n+1> = U_n|psi_n>
    psi_reshaped = psi.reshape(Nx, Ny)

    # Extraire la fonction d'onde sur l'écran
    for i in range(len(screen_index_x)) :
        psi_at_screen = psi_reshaped[:, screen_index_x[i]]  # On prend la colonne pour x = position de lecran
        profiles_on_screen[i].append(np.abs(psi_at_screen)**2) # Ajouter le profil sur l'écran à la liste

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
    E_inst = compute_energy(psi, H)
    E_moy += E_inst
    # mesure des moyennes de norme
    normalisation_inst = np.sum(np.abs(psi)**2)
    normalisation_moy += normalisation_inst
    if (n % mesure_interval == 0) :
        # mesure de l'énergie
        E_moy /= mesure_interval # fait une moyenne
        energy_evolution.append(E_moy)
        E_moy = 0
        # mesure de la norme
        normalisation_moy /= mesure_interval
        normalisation_evolution.append(normalisation_moy)
        normalisation_moy = 0

        # temps ou on a fait la mesure
        mesure_time.append(n*dt)
        
    if n % animation_interval  == 0 : #on sauvegarde les positions pour les annimées tout les "animation_interval" pas 
        probability_for_animation.append(np.abs(psi_reshaped)**2)
        print(f"image {n/animation_interval}/{Nt/animation_interval}")

print("energie finale = ", energy_evolution[-1], "norme finale = ", normalisation_evolution[-1])
print("Simulation terminée")

# lance l'animation
if save_frames and len(probability_for_animation) > 1 :
    animate_wavefunction_2D(probability_for_animation, Vxy)

# temps que met le programme
end_time = time.perf_counter()
execution_time = end_time - start_time
print(f"Programme exécuté en : {convert(execution_time)}")

# plot le potentiel en 3D
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Vxy, cmap='Blues')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('V(x, y)')
ax.set_title('Potentiel en 2D')
plt.show()

# Affichage des bandes d'interférence
for i in range(len(SCREEN_POSITION)):
    final_profile = np.array(profiles_on_screen[i])  # Tableau de forme (Nt, Ny)

    # Création de la figure avec 3 sous-graphiques (1 colonne, 3 lignes)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5)) 

    # 1. Affichage de l'évolution temporelle des bandes d'interférence
    im1 = axes[0].imshow(final_profile, cmap="inferno", aspect="auto", extent=[-Ly/2, Ly/2, 0, Nt * dt])
    axes[0].set_xlabel("y")
    axes[0].set_ylabel("Temps")
    axes[0].set_title(f"Évolution des bandes d'interférence sur l'écran en x = {SCREEN_POSITION[i]}")
    fig.colorbar(im1, ax=axes[0], label="|ψ|²")  # Ajouter une barre de couleur

    # 2. Tracé de la somme des probabilités sur l'écran
    axes[1].plot(np.linspace(-Ly/2, Ly/2, Nx), np.sum(final_profile, axis=0), color='r')
    axes[1].set_title(f"Bandes d'interférence sur l'écran en x = {SCREEN_POSITION[i]}")
    axes[1].set_xlabel("Position en y")
    axes[1].set_ylabel("Somme des probabilités au cours du temps")
    axes[1].grid(True)

    # 3. Projection 2D de la somme des probabilités
    Nz = Ny
    dz = dy
    screen = np.tile(np.sum(final_profile, axis=0), (Nz, 1))  # Réplication correcte

    im3 = axes[2].imshow(screen, cmap="inferno", aspect="auto", extent=[-Ly/2, Ly/2, 0, Nz * dz])
    axes[2].set_xlabel("y")
    axes[2].set_ylabel("z")
    axes[2].set_title(f"L'écran en x = {SCREEN_POSITION[i]}")
    fig.colorbar(im3, ax=axes[2], label="Intensité")

    plt.tight_layout()  # Ajuste automatiquement les espacements
    plt.show()

if not(ABSORPTION) : #car si on a l absorption l energie et la norme n est pas constante or avec ce graphique on veut montrer que la méthode est stable
    # Tracer l'évolution de l'énergie
    plt.plot((np.arange(0, Nt, mesure_interval) * dt)[1:], energy_evolution[1:]) #on retire la premiere valeur qui est a 0 a cause de son initialisation
    plt.xlabel("Temps (s)")
    plt.ylabel("Énergie")
    plt.title("Évolution de l'énergie au cours de la simulation")
    plt.show()

    # Tracer l'évolution de la norme
    plt.plot((np.arange(0, Nt, mesure_interval) * dt)[1:], normalisation_evolution[1:])
    plt.xlabel("Temps (s)")
    plt.ylabel("|psi|²")
    plt.title("Évolution de la norme de psi au cours de la simulation")
    plt.show()
