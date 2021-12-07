from random import randint
from math import cos, sin, radians
import time

#######        Variables globales        #######

#Genetiques
NB_POP = 200
NB_GENES = 80
ELITE = 20          #le top 30% est sélectionné
RANDOM = 0          #20% sont choisi aléatoirement pour rester
MUTATION = 5        #Sur la nouvelle population, 6% des genes totaux sont modifiées
numéro_génération = 0
#Affichage
DT = 2
#Sol de Mars
x_init = 2500
y_init = 2700
l = [(0, 100), (1000, 500), (1500, 1500), (3000, 1000), (4000, 150), (5500, 150), (6999, 800)]


#######        Fonctions utiles        #######
def generation_genes():
    """Genere 100 genes qui modélise une puissance et une orientation aléatoire"""
    return [(randint(-1, 1), randint(-15, 15)) for i in range (NB_GENES)]

def landing_area(l):
    """On calcule les coordonnées de la zone d'atterissage"""
    for i in range (len(l) - 1):
        if l[i][1] == l[i+1][1] and (l[i+1][0] - l[i][0]) > 1000:
            return l[i][0], l[i+1][0], l[i][1]

#Global coordinates of the landing area
x1_land, x2_land, y_land = landing_area(l)

def in_landing_area(x, y):
    """Verifie si le module est en zone d'atterissage ou non"""
    if x >= x1_land and x <= x2_land and y > y_land:
        return True

def calcul_score(d):
    """Calcul le score en fonction de la distance de la zone d'atterissage
    si la distance est nulle, le score est de 52
    plus la distance est grande plus le score sera faible
    Au dessus de 3500 le score est nulle"""
    if d > 3500**2:
        return 0
    else:
        return 28 - (d / (3500**2)) * 28

def score(chrom):
    """Score de cet échantillon en fonction:
    Le score maximal est de 100 
    52 points pour la distance de la zone d'atterissage
    Si il est dans la zone d'atterissage:
    16 si son orientation est correct
    16 si son vecteur vitesse suivant x est suffissament faible
    16 si son vecteur vitesse suivant y est suffissament faible"""
    if chrom.x < x1_land:
        return calcul_score(abs((x1_land - chrom.x)**2 + (y_land - chrom.y)**2))
    elif chrom.x > x2_land:
        return calcul_score(abs((x2_land - chrom.x)**2 + (y_land - chrom.y)**2))
    else:
        s = 28
        #Vitesse verticale
        if abs(chrom.v_y) <= 40: 
            s += 24
            #Vitesse horizontale
            if abs(chrom.v_x) <= 20:
                s +=  24
                #Inclinaison
                s += (24 - abs(chrom.inclinaison) * 0.27)
            elif abs(chrom.v_x) <= 80:
                s += 24 - (abs(chrom.v_x) - 20) * 0.4
        elif abs(chrom.v_y) <= 100:
            s += 24 - (abs(chrom.v_y) - 40) * 0.4
    return s       

def collision(x, y, x_pre, y_pre):
    """"Calcul de la distance avec le sol pour savoir si le module à atteint le sol"""
    #On calcul premièrement les intervalles pour le segment de la trajectoire
    a = min(x_pre, x)
    b = max(x_pre, x)
    c = min(y_pre, y)
    d = max(y_pre, y)
    #Pour chaque segment du sol on calcul pour savoir si il y a une intersection ou non
    for i in range (len(l) - 1):
        x1, y1 = l[i]
        x2, y2 = l[i+1]

        e = min(x1, x2)
        f = max(x1, x2)
        g = min(y1, y2)
        h = max(y1, y2)

        if b < e or a > f or c > h or d < g:
            continue
        else:
            if x == x_pre:
                if x1 == x2:
                    if x1 != x:
                        continue
                    else:
                        return True
                else:
                    t = (x - e) * (y2 - y1) / (x2 - x1)
                    if t <= d and t >= c:
                        return True
                    else:
                        continue
            else:
                if x1 == x2:
                    t = (x1 - a) * (y - y_pre) / (x - x_pre)
                    if t <= h and t >= g:
                        return True
                    else:
                        continue
                else:
                    a1 = (y2 - y1) / (x2 - x1)
                    a2 = (y - y_pre) / (x - x_pre)
                    b1 = y1 - a1 * x1
                    b2 = y - a2 * x
                    if a1 == a2:
                        continue
                    else:
                        xa = (b2 - b1) / (a1 - a2)
                        if xa < max(a, e) or xa > min(b, f):
                            continue
                        else:
                            return True
    return False
                   
def arret(chrom):
    """Verifie la distance du module par rapport au sol
    si il est trop proche, il s'est arrêté plus ou moins correctement"""

    if chrom.x > 7000 or chrom.x < 0 or chrom.y > 3000 or chrom.y < 0:
        return True, score(chrom)
    else:
        if collision(chrom.x, chrom.y, chrom.previous_x, chrom.previous_y):
            return True, score(chrom)
        else:
            return False, 0


#######        Classes        #######
class Chromo:
    """Represente chaque chromosome soumis à l'attraction de Mars et à la force
    de ses réacteurs"""
    def __init__(self, genes = None):
        self.previous_x = x_init
        self.previous_y = y_init
        self.first = True
        self.x = x_init
        self.y = y_init
        self.v_x = 0
        self.v_y = 0
        self.a_x = 0
        self.a_y = 0
        self.iteration = 0
        self.power = 0
        self.inclinaison = 0
        self.alive = True
        self.genes=generation_genes() if genes==None else genes

    def update(self):
        a = self.genes[self.iteration]
        self.power += a[0]
        self.inclinaison += a[1]
        if self.power<0: self.power = 0
        elif self.power>4: self.power = 4

        if self.inclinaison<-90: self.inclinaison = -90
        elif self.inclinaison>90: self.inclinaison = 90

        self.iteration += 1

        self.a_x = - self.power * sin(radians(self.inclinaison))
        self.a_y = self.power * cos(radians(self.inclinaison)) - 3.711
        self.v_x += DT * self.a_x
        self.v_y += DT * self.a_y
        if self.first:
            self.first = False
        else:
            self.previous_x = self.x
            self.previous_y = self.y
        self.x += 0.5 * DT * DT * self.a_x + DT * self.v_x
        self.y += 0.5 * DT * DT * self.a_y + DT * self.v_y

    def reset(self):
        self.previous_x = x_init
        self.previous_y = y_init
        self.x = x_init
        self.y = y_init
        self.v_x = 0
        self.v_y = 0
        self.iteration = 0
        self.power = 0
        self.inclinaison = 0
        self.alive = True
        self.first = True


class Population:
    def __init__(self):
        self.pop = []
        self.resultat = {}
        self.count_dead = 0
    
    def creation(self):
        self.pop = [Chromo() for i in range (NB_POP)]

    def update(self):
        for i in range (NB_POP):
            chrom = self.pop[i]
            if chrom.alive:
                b, s = arret(chrom)
                if not b:
                    chrom.update()
                else:
                    self.resultat[i] = s
                    chrom.alive = False
                    self.count_dead += 1

    def vainqueur(self):
        values = list(self.resultat.values())
        if 100 in values:
            return True, self.pop[values.index(100)].genes
        else:
            return False, 0

    def nouvelle_génération(self):
        #Trier le dico pour en sortir les meilleurs éléments
        result = sorted(self.resultat.items(), key=lambda x:x[1], reverse=True)
        #30% des meilleurs des meilleurs des meilleurs
        pourcentage_elite = int(NB_POP * ELITE / 100)
        nouvelle_gen_index = [a[0] for a in result[:pourcentage_elite]]
        reste = [a[0] for a in result[pourcentage_elite:]]
        #20% de rattrapage
        while len(nouvelle_gen_index) < int(NB_POP * (ELITE + RANDOM) / 100):
            i = randint(0, len(reste) - 1)
            a = reste[i]
            del reste[i]
            nouvelle_gen_index.append(a)
        #Récupération des chromosomes avec les index choisis
        nouvelle_gen = []
        for i in nouvelle_gen_index:
            chrom = self.pop[i]
            chrom.reset()
            nouvelle_gen.append(chrom)
        #Reproduction
        while len(nouvelle_gen) < NB_POP:
            i = randint(0, len(nouvelle_gen) - 1)
            j = randint(0, len(nouvelle_gen) - 1)
            new_genes = nouvelle_gen[i].genes[:(NB_GENES//2)] + nouvelle_gen[j].genes[(NB_GENES//2):]
            nouvelle_gen.append(Chromo(new_genes))
        #Mutation
        for i in range (int(NB_POP * NB_GENES * MUTATION / 100)):
            j = randint(0, NB_POP - 1)
            k = randint(0, NB_GENES - 1)
            nouvelle_gen[j].genes[k] = (randint(-1, 1), randint(-15, 15))
        #Le grand remplacement
        self.pop = nouvelle_gen
        self.resultat = {}
        self.count_dead = 0


def main():
    début = time.time()
    global numéro_génération

    #Creer la population
    populace = Population()
    populace.creation()

    #Affichage de la fenetre creee
    while not populace.vainqueur()[0]:
        while populace.count_dead != NB_POP:
            populace.update()
            
        a, b = populace.vainqueur()
        if a:
            print(b)
            print(numéro_génération)
            print(time.time() - début)
            break
        else:
            populace.nouvelle_génération()
            numéro_génération += 1

if __name__ == "__main__":
    main()